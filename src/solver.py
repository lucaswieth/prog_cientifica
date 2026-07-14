import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import logging

from malha import Malha
from advection_schemes import AdvectionScheme
from boundary_conditions import BoundaryConditions
from config import SimulationConfig

class TransportSolver:
    """
    Solves the 2D Advection-Diffusion-Reaction equation for a scalar phi 
    on a structured grid using Finite Volume Method (FVM) and Implicit Euler.
    """

    def __init__(self, grid: Malha, scheme: AdvectionScheme, config: SimulationConfig, logger: logging.Logger):
        self.grid = grid
        self.scheme = scheme
        self.config = config
        self.logger = logger
        self._peclet_warning_emitted = False

    def _get_index(self, i: int, j: int) -> int:
        """ Maps 2D coordinates to a 1D global index. """
        return i * self.grid.Ny + j

    def _check_stability(self, F: float, D: float) -> None:
        """
        Checks the local Peclet number for stability and emits a warning if needed.
        Pe = |F| / D
        """
        if D > 0.0:
            peclet = abs(F) / D
            if peclet > self.config.critical_peclet and not self._peclet_warning_emitted:
                self.logger.warning(
                    f"Local Peclet number ({peclet:.2f}) exceeds critical value "
                    f"({self.config.critical_peclet}). Central Difference Scheme might produce non-physical oscillations."
                )
                self._peclet_warning_emitted = True

    def assemble_system(self, bc: BoundaryConditions, k: float = 0.0, scalar_type: str = None, 
                        phi2_field: np.ndarray = None, dt: float = None, phi_old: np.ndarray = None) -> tuple[lil_matrix, np.ndarray]:
        """
        Assembles the global linear system A * phi = b.
        Uses Implicit Euler if dt and phi_old are provided.
        """
        Nx = self.grid.Nx
        Ny = self.grid.Ny
        N_total = Nx * Ny

        A = lil_matrix((N_total, N_total))
        b = np.zeros(N_total)

        dx = 1.0 / Nx
        dy = 1.0 / Ny
        V = dx * dy

        # Diffusion conductances
        D_e = self.config.gamma * dy / dx
        D_w = self.config.gamma * dy / dx
        D_n = self.config.gamma * dx / dy
        D_s = self.config.gamma * dx / dy

        for i in range(Nx):
            for j in range(Ny):
                P = self._get_index(i, j)

                # Convective fluxes
                F_e = self.grid.u_faces_x[i+1, j] * dy
                F_w = self.grid.u_faces_x[i, j] * dy
                F_n = self.grid.v_faces_y[i, j+1] * dx
                F_s = self.grid.v_faces_y[i, j] * dx

                # Check stability for all faces
                self._check_stability(F_e, D_e)
                self._check_stability(F_w, D_w)
                self._check_stability(F_n, D_n)
                self._check_stability(F_s, D_s)

                a_P = 0.0

                # --- EAST FACE ---
                if i < Nx - 1:
                    E = self._get_index(i+1, j)
                    a_E, a_p_contrib = self.scheme.calculate_coefficients(D_e, F_e, is_positive_direction=True)
                    A[P, E] = -a_E
                    a_P += a_p_contrib
                else:
                    a_P += 2.0 * D_e
                    b[P] += (2.0 * D_e - F_e) * bc.east

                # --- WEST FACE ---
                if i > 0:
                    W = self._get_index(i-1, j)
                    a_W, a_p_contrib = self.scheme.calculate_coefficients(D_w, F_w, is_positive_direction=False)
                    A[P, W] = -a_W
                    a_P += a_p_contrib
                else:
                    a_P += 2.0 * D_w
                    b[P] += (2.0 * D_w + F_w) * bc.west

                # --- NORTH FACE ---
                if j < Ny - 1:
                    N = self._get_index(i, j+1)
                    a_N, a_p_contrib = self.scheme.calculate_coefficients(D_n, F_n, is_positive_direction=True)
                    A[P, N] = -a_N
                    a_P += a_p_contrib
                else:
                    a_P += 2.0 * D_n
                    b[P] += (2.0 * D_n - F_n) * bc.north

                # --- SOUTH FACE ---
                if j > 0:
                    S = self._get_index(i, j-1)
                    a_S, a_p_contrib = self.scheme.calculate_coefficients(D_s, F_s, is_positive_direction=False)
                    A[P, S] = -a_S
                    a_P += a_p_contrib
                else:
                    a_P += 2.0 * D_s
                    b[P] += (2.0 * D_s + F_s) * bc.south

                # --- TRANSIENT TERM (Implicit Euler) ---
                if dt is not None and phi_old is not None:
                    a_P += V / dt
                    b[P] += (V / dt) * phi_old[i, j]

                # --- REACTION TERM (Source/Sink) ---
                if k > 0.0:
                    if scalar_type == 'phi2':
                        # phi2 is consumed: S_P = -k * phi2 * V => goes to LHS as +k*V
                        A[P, P] += k * V
                    elif scalar_type == 'phi1':
                        # phi1 is produced from phi2: S_C = k * phi2 * V => goes to RHS
                        if phi2_field is not None:
                            b[P] += k * phi2_field[i, j] * V
                        else:
                            raise ValueError("phi2_field is required to solve phi1 with k > 0.")

                # Assign central coefficient
                A[P, P] += a_P

        return A, b

    def solve_transient(self, bc: BoundaryConditions, t_final: float, dt: float, 
                        k: float = 0.0, scalar_type: str = None, phi2_hist: list[np.ndarray] = None) -> list[np.ndarray]:
        """
        Solves the transient problem from t=0 (phi=0) to t_final.
        Returns a history of the 2D scalar fields over time.
        """
        tempos = np.arange(0, t_final + dt, dt)
        historico = []
        
        phi_atual = np.zeros((self.grid.Nx, self.grid.Ny))
        historico.append(phi_atual.copy())
        
        for n, t in enumerate(tempos[1:], start=1):
            phi2_field = phi2_hist[n] if (phi2_hist is not None and len(phi2_hist) > n) else None
            
            A_lil, b = self.assemble_system(bc, k=k, scalar_type=scalar_type, phi2_field=phi2_field, dt=dt, phi_old=phi_atual)
            A_csr = A_lil.tocsr()
            phi_1d = spsolve(A_csr, b)
            
            phi_atual = phi_1d.reshape((self.grid.Nx, self.grid.Ny))
            historico.append(phi_atual.copy())
            
        return historico
