import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class PostProcessor:
    """
    Handles plotting and animation of the CFD results.
    """

    @staticmethod
    def plotar_resultados(phi_2d: np.ndarray, grid, titulo: str, save_path: str = None) -> None:
        """ 
        Renders the final contours with overlaid streamlines of the grid velocities.
        """
        X, Y = np.meshgrid(grid.x_centros, grid.y_centros, indexing='ij')

        fig, ax = plt.subplots(figsize=(7.5, 6))
        fixed_levels = np.linspace(-0.01, 1.01, 25)
        contour = ax.contourf(X, Y, phi_2d, levels=fixed_levels, cmap='plasma', vmin=-0.01, vmax=1.01)
        cbar = fig.colorbar(contour, ax=ax)
        cbar.set_label('Valor de $\\phi$', fontsize=12)

        # Reconstruct velocities at cell centers for streamlines
        U_centros = 0.5 * (grid.u_faces_x[:-1, :] + grid.u_faces_x[1:, :])
        V_centros = 0.5 * (grid.v_faces_y[:, :-1] + grid.v_faces_y[:, 1:])
        
        ax.streamplot(grid.x_centros, grid.y_centros, U_centros.T, V_centros.T, color='white', linewidth=0.8, density=1.0, arrowsize=1.0)

        ax.set_title(titulo, fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel('Coordenada X (m)', fontsize=12)
        ax.set_ylabel('Coordenada Y (m)', fontsize=12)
        ax.set_aspect('equal')
        ax.grid(True, linestyle=':', alpha=0.5)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')

        plt.close()

    @staticmethod
    def gerar_animacao_temporal(historico: list[np.ndarray], grid, dt: float, titulo: str, save_path: str) -> None:
        """ 
        Generates and saves a time-series animation (.mp4 or .gif) from transient history.
        """
        fig, ax = plt.subplots(figsize=(7.5, 6))
        X, Y = np.meshgrid(grid.x_centros, grid.y_centros, indexing='ij')
        
        vmin, vmax = -0.01, 1.01
        
        fixed_levels = np.linspace(vmin, vmax, 25)
            
        U_centros = 0.5 * (grid.u_faces_x[:-1, :] + grid.u_faces_x[1:, :])
        V_centros = 0.5 * (grid.v_faces_y[:, :-1] + grid.v_faces_y[:, 1:])
        
        # Cria a colorbar usando fixed_levels para evitar que um frame zerado corrompa a escala
        contour_init = ax.contourf(X, Y, historico[-1], levels=fixed_levels, cmap='plasma', vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(contour_init, ax=ax)
        cbar.set_label('Valor de $\\phi$', fontsize=12)
        
        def update(frame):
            ax.clear()
            ax.contourf(X, Y, historico[frame], levels=fixed_levels, cmap='plasma', vmin=vmin, vmax=vmax)
            ax.streamplot(grid.x_centros, grid.y_centros, U_centros.T, V_centros.T, color='white', linewidth=0.5, density=1.0)
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
            raise ValueError("File format must be .mp4 or .gif")
            
        plt.close()

    @staticmethod
    def plotar_corte_1d(phi_uds: np.ndarray, phi_cds: np.ndarray, grid, title: str, save_path: str = None) -> None:
        """
        Plots a 1D vertical cut (fixed x in the middle) comparing UDS and CDS.
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
            
        plt.close()
