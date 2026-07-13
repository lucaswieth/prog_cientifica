import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class TransportSolver:
    """
    Classe responsável por resolver a equação de advecção-difusão-reação bidimensional (2D)
    para o transporte de um escalar phi em uma malha estruturada. Suporta formulações
    estacionárias e transientes (Euler Implícito).
    """

    def __init__(self, grid, scheme: str = 'UDS', gamma: float = 0.001):
        self.grid = grid
        self.scheme = scheme.upper()
        self.gamma = gamma

        if self.scheme not in ['UDS', 'CDS']:
            raise ValueError("O esquema advectivo deve ser 'UDS' ou 'CDS'")

    def _get_index(self, i: int, j: int) -> int:
        return i * self.grid.Ny + j

    def montar_sistema(self, phi_bc: dict, k: float = 0.0, scalar_type: str = None, 
                       phi2_field: np.ndarray = None, dt: float = None, phi_old: np.ndarray = None):
        """
        Monta o sistema linear global A * phi = b discretizado por Volumes Finitos.
        Se dt e phi_old forem fornecidos, ativa o termo transiente (Euler Implícito).
        """
        Nx = self.grid.Nx
        Ny = self.grid.Ny
        N_total = Nx * Ny

        A = lil_matrix((N_total, N_total))
        b = np.zeros(N_total)

        dx = 1.0 / Nx
        dy = 1.0 / Ny
        V = dx * dy

        D_e = self.gamma * dy / dx
        D_w = self.gamma * dy / dx
        D_n = self.gamma * dx / dy
        D_s = self.gamma * dx / dy

        for i in range(Nx):
            for j in range(Ny):
                P = self._get_index(i, j)

                F_e = self.grid.u_faces[i+1, j] * dy
                F_w = self.grid.u_faces[i, j] * dy
                F_n = self.grid.v_faces[i, j+1] * dx
                F_s = self.grid.v_faces[i, j] * dx

                a_E = 0.0
                a_W = 0.0
                a_N = 0.0
                a_S = 0.0
                a_P = 0.0

                # FACE LESTE
                if i < Nx - 1:
                    E = self._get_index(i+1, j)
                    if self.scheme == 'UDS':
                        a_E = D_e + max(-F_e, 0.0)
                    else:
                        a_E = D_e - 0.5 * F_e
                    A[P, E] = -a_E
                    a_P += a_E + F_e
                else:
                    a_P += 2.0 * D_e
                    b[P] += (2.0 * D_e - F_e) * phi_bc['east']

                # FACE OESTE
                if i > 0:
                    W = self._get_index(i-1, j)
                    if self.scheme == 'UDS':
                        a_W = D_w + max(F_w, 0.0)
                    else:
                        a_W = D_w + 0.5 * F_w
                    A[P, W] = -a_W
                    a_P += a_W - F_w
                else:
                    a_P += 2.0 * D_w
                    b[P] += (2.0 * D_w + F_w) * phi_bc['west']

                # FACE NORTE
                if j < Ny - 1:
                    N = self._get_index(i, j+1)
                    if self.scheme == 'UDS':
                        a_N = D_n + max(-F_n, 0.0)
                    else:
                        a_N = D_n - 0.5 * F_n
                    A[P, N] = -a_N
                    a_P += a_N + F_n
                else:
                    a_P += 2.0 * D_n
                    b[P] += (2.0 * D_n - F_n) * phi_bc['north']

                # FACE SUL
                if j > 0:
                    S = self._get_index(i, j-1)
                    if self.scheme == 'UDS':
                        a_S = D_s + max(F_s, 0.0)
                    else:
                        a_S = D_s + 0.5 * F_s
                    A[P, S] = -a_S
                    a_P += a_S - F_s
                else:
                    a_P += 2.0 * D_s
                    b[P] += (2.0 * D_s + F_s) * phi_bc['south']

                # --- Termo Transiente (Euler Implícito) ---
                if dt is not None and phi_old is not None:
                    a_P += V / dt
                    b[P] += (V / dt) * phi_old[i, j]

                A[P, P] = a_P

                # --- Termo Fonte de Reação (k) ---
                if k > 0.0:
                    if scalar_type == 'phi2':
                        A[P, P] += k * V
                    elif scalar_type == 'phi1':
                        if phi2_field is not None:
                            b[P] += k * phi2_field[i, j] * V
                        else:
                            raise ValueError("Campo phi2_field é obrigatório para resolver phi1 com k > 0.")

        return A, b

    def resolver(self, phi_bc: dict, k: float = 0.0, scalar_type: str = None, phi2_field: np.ndarray = None) -> np.ndarray:
        """ Resolve o problema de forma estacionária pura. """
        A_lil, b = self.montar_sistema(phi_bc, k, scalar_type, phi2_field)
        A_csr = A_lil.tocsr()
        phi_1d = spsolve(A_csr, b)
        phi_2d = phi_1d.reshape((self.grid.Nx, self.grid.Ny))
        return phi_2d

    def resolver_transiente(self, phi_bc: dict, t_final: float, dt: float, k: float = 0.0, scalar_type: str = None, phi2_hist: list = None) -> list:
        """
        Resolve a equação no tempo usando Euler Implícito a partir da condição inicial phi=0.
        Retorna uma lista de campos 2D correspondentes a cada passo de tempo (incluindo t=0).
        """
        tempos = np.arange(0, t_final + dt, dt)
        historico = []
        
        phi_atual = np.zeros((self.grid.Nx, self.grid.Ny))
        historico.append(phi_atual.copy())
        
        for n, t in enumerate(tempos[1:], start=1):
            # Para a reação, phi1 precisa da concentração de phi2 no mesmo instante de tempo
            phi2_field = phi2_hist[n] if (phi2_hist is not None and len(phi2_hist) > n) else None
            
            A_lil, b = self.montar_sistema(phi_bc, k=k, scalar_type=scalar_type, phi2_field=phi2_field, dt=dt, phi_old=phi_atual)
            A_csr = A_lil.tocsr()
            phi_1d = spsolve(A_csr, b)
            phi_atual = phi_1d.reshape((self.grid.Nx, self.grid.Ny))
            historico.append(phi_atual.copy())
            
        return historico

    def plotar_resultados(self, phi_2d: np.ndarray, titulo: str, save_path: str = None):
        """ Renderiza os contornos finais com as streamlines da malha sobrepostas. """
        X, Y = np.meshgrid(self.grid.x_centros, self.grid.y_centros, indexing='ij')

        fig, ax = plt.subplots(figsize=(7.5, 6))
        contour = ax.contourf(X, Y, phi_2d, levels=25, cmap='plasma')
        cbar = fig.colorbar(contour, ax=ax)
        cbar.set_label('Valor de $\\phi$', fontsize=12)

        # Streamlines sobrepostas
        U_centros = 0.5 * (self.grid.u_faces[:-1, :] + self.grid.u_faces[1:, :])
        V_centros = 0.5 * (self.grid.v_faces[:, :-1] + self.grid.v_faces[:, 1:])
        ax.streamplot(self.grid.x_centros, self.grid.y_centros, U_centros.T, V_centros.T, color='white', linewidth=0.8, density=1.0, arrowsize=1.0)

        ax.set_title(titulo, fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel('Coordenada X (m)', fontsize=12)
        ax.set_ylabel('Coordenada Y (m)', fontsize=12)
        ax.set_aspect('equal')
        ax.grid(True, linestyle=':', alpha=0.5)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Gráfico salvo com sucesso em: {save_path}")

        plt.close()

    def gerar_animacao_temporal(self, historico: list, dt: float, titulo: str, save_path: str):
        """ Gera e salva animação em .mp4 ou .gif a partir do histórico de campos transientes. """
        fig, ax = plt.subplots(figsize=(7.5, 6))
        X, Y = np.meshgrid(self.grid.x_centros, self.grid.y_centros, indexing='ij')
        
        vmax = max([np.max(phi) for phi in historico])
        vmin = min([np.min(phi) for phi in historico])
        if vmax == vmin:
            vmax = vmin + 1.0 
            
        U_centros = 0.5 * (self.grid.u_faces[:-1, :] + self.grid.u_faces[1:, :])
        V_centros = 0.5 * (self.grid.v_faces[:, :-1] + self.grid.v_faces[:, 1:])
        
        def update(frame):
            ax.clear()
            ax.contourf(X, Y, historico[frame], levels=25, cmap='plasma', vmin=vmin, vmax=vmax)
            ax.streamplot(self.grid.x_centros, self.grid.y_centros, U_centros.T, V_centros.T, color='white', linewidth=0.5, density=1.0)
            ax.set_title(f"{titulo} (t = {frame*dt:.2f} s)", fontsize=14, fontweight='bold', pad=15)
            ax.set_xlabel('Coordenada X (m)', fontsize=12)
            ax.set_ylabel('Coordenada Y (m)', fontsize=12)
            ax.set_aspect('equal')
            ax.grid(True, linestyle=':', alpha=0.5)

        anim = animation.FuncAnimation(fig, update, frames=len(historico), interval=100)
        
        if save_path.endswith('.mp4'):
            anim.save(save_path, writer='ffmpeg', dpi=100)
        elif save_path.endswith('.gif'):
            anim.save(save_path, writer='pillow', dpi=100)
        else:
            raise ValueError("O formato do arquivo deve ser .mp4 ou .gif")
            
        plt.close()

    @staticmethod
    def plotar_corte_1d(phi_uds: np.ndarray, phi_cds: np.ndarray, grid, title: str, save_path: str = None):
        """
        Plota um corte vertical (x fixo no meio) das concentrações, 
        comparando esquemas UDS e CDS.
        """
        fig, ax = plt.subplots(figsize=(6, 5))
        
        i_meio = grid.Nx // 2
        y_coords = grid.y_centros
        
        ax.plot(y_coords, phi_uds[i_meio, :], 'b-o', label='UDS', markersize=4)
        ax.plot(y_coords, phi_cds[i_meio, :], 'r-s', label='CDS', markersize=4)
        
        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
        ax.set_xlabel('Coordenada Y (m)')
        ax.set_ylabel('Concentração $\\phi$ em X = 0.5')
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend()
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            
        # Não chamamos plt.show() para evitar travamentos
        plt.close()
