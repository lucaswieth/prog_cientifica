import numpy as np
from scipy.interpolate import griddata

class Malha:
    """
    Finite Volume Method (FVM) 2D Grid Generator and Interpolator.
    Manages geometric cell centers, faces, and interpolates an external frozen velocity field onto them.
    """

    def __init__(self, x_orig: np.ndarray, y_orig: np.ndarray, u_orig: np.ndarray, v_orig: np.ndarray, 
                 Nx: int, Ny: int, interp_method: str = 'linear'):
        """
        Initializes the grid and interpolates the given velocity field.

        Args:
            x_orig (np.ndarray): Original x coordinates of the dataset.
            y_orig (np.ndarray): Original y coordinates of the dataset.
            u_orig (np.ndarray): Original x-velocity of the dataset.
            v_orig (np.ndarray): Original y-velocity of the dataset.
            Nx (int): Number of control volumes in x.
            Ny (int): Number of control volumes in y.
            interp_method (str): SciPy interpolation method ('linear', 'cubic', 'nearest').
        """
        self.Nx = Nx
        self.Ny = Ny
        self.interp_method = interp_method
        
        self._points_orig = np.column_stack((x_orig, y_orig))
        self._u_orig = u_orig
        self._v_orig = v_orig

        self._build_grid()
        self._interpolate_velocities()

    def _build_grid(self) -> None:
        """
        Generates face coordinates and cell centers for a domain [0, 1] x [0, 1].
        """
        self.x_faces = np.linspace(0.0, 1.0, self.Nx + 1)
        self.y_faces = np.linspace(0.0, 1.0, self.Ny + 1)

        self.x_centros = 0.5 * (self.x_faces[:-1] + self.x_faces[1:])
        self.y_centros = 0.5 * (self.y_faces[:-1] + self.y_faces[1:])

        # Meshgrids for 2D interpolation
        self.X_face_x, self.Y_face_x = np.meshgrid(self.x_faces, self.y_centros, indexing='ij')
        self.X_face_y, self.Y_face_y = np.meshgrid(self.x_centros, self.y_faces, indexing='ij')
        self.X_c, self.Y_c = np.meshgrid(self.x_centros, self.y_centros, indexing='ij')

    def _interp_2d_with_fallback(self, values_src: np.ndarray, x_query: np.ndarray, y_query: np.ndarray) -> np.ndarray:
        """
        Interpolates a scalar field with a nearest-neighbor fallback for NaNs at boundaries.
        """
        grid_val = griddata(self._points_orig, values_src, (x_query, y_query), method=self.interp_method)

        nan_mask = np.isnan(grid_val)
        if np.any(nan_mask):
            grid_val_nearest = griddata(self._points_orig, values_src, (x_query, y_query), method='nearest')
            grid_val[nan_mask] = grid_val_nearest[nan_mask]

        return grid_val

    def _interpolate_velocities(self) -> None:
        """
        Interpolates the original velocity field onto the geometric grid points.
        """
        # East/West interfaces
        self.u_faces_x = self._interp_2d_with_fallback(self._u_orig, self.X_face_x, self.Y_face_x)
        self.v_faces_x = self._interp_2d_with_fallback(self._v_orig, self.X_face_x, self.Y_face_x)

        # North/South interfaces
        self.u_faces_y = self._interp_2d_with_fallback(self._u_orig, self.X_face_y, self.Y_face_y)
        self.v_faces_y = self._interp_2d_with_fallback(self._v_orig, self.X_face_y, self.Y_face_y)

        # Cell centers
        self.u_centros = self._interp_2d_with_fallback(self._u_orig, self.X_c, self.Y_c)
        self.v_centros = self._interp_2d_with_fallback(self._v_orig, self.X_c, self.Y_c)
