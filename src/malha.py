import os
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

class Malha:
    """
    Classe responsável pela geometria da malha bidimensional centrada na célula (FVM 2D)
    e pelo gerenciamento do campo de velocidades (frozen flow) por interpolação de dados.
    """

    def __init__(self, data_filepath: str, Nx: int, Ny: int, interp_method: str = 'linear'):
        """
        Inicializa a malha física [0, 1] x [0, 1] e interpola os dados de escoamento.

        Parameters
        ----------
        data_filepath : str
            Caminho do arquivo LidDrivenCavityFlow_Re_1000.dat.
        Nx : int
            Número de volumes de controle na direção x.
        Ny : int
            Número de volumes de controle na direção y.
        interp_method : str, opcional
            Método de interpolação do SciPy ('linear', 'cubic', 'nearest'). O padrão é 'linear'.
        """
        self.Nx = Nx
        self.Ny = Ny
        self.data_filepath = data_filepath
        self.interp_method = interp_method

        # 1. Leitura e limpeza dos dados
        self._load_and_clean_data()

        # 2. Construção da malha geométrica (FVM Cell-Centered)
        self._build_grid()

        # 3. Interpolação de velocidades
        self._interpolate_velocities()

    def _load_and_clean_data(self):
        """
        Lê os dados originais de escoamento e converte-os para arrays do NumPy.
        """
        if not os.path.exists(self.data_filepath):
            raise FileNotFoundError(f"Arquivo de dados não encontrado: {self.data_filepath}")

        # Leitura via pandas
        df = pd.read_csv(self.data_filepath)

        # Remoção de espaços extras nos cabeçalhos
        df.columns = df.columns.str.strip()

        # Armazenamento e conversão direta para arrays do NumPy
        self._x_orig = df['x-coordinate'].to_numpy()
        self._y_orig = df['y-coordinate'].to_numpy()
        self._u_orig = df['x-velocity'].to_numpy()
        self._v_orig = df['y-velocity'].to_numpy()

        # Matriz de pontos de origem (N_pontos, 2)
        self._points_orig = np.column_stack((self._x_orig, self._y_orig))

    def _build_grid(self):
        """
        Gera as coordenadas das faces e centros dos volumes de controle no domínio [0, 1] x [0, 1].
        """
        # Coordenadas das faces (arrays 1D)
        self._x_faces = np.linspace(0.0, 1.0, self.Nx + 1)
        self._y_faces = np.linspace(0.0, 1.0, self.Ny + 1)

        # Coordenadas dos centros das células (média das faces adjacentes)
        self._x_centros = 0.5 * (self._x_faces[:-1] + self._x_faces[1:])
        self._y_centros = 0.5 * (self._y_faces[:-1] + self._y_faces[1:])

        # Meshgrids úteis para a interpolação bi-dimensional
        # 1. Faces normais a x (interfaces de fluxo horizontal - East/West): coordenadas (x_faces, y_centros)
        self.X_face_x, self.Y_face_x = np.meshgrid(self._x_faces, self._y_centros, indexing='ij')

        # 2. Faces normais a y (interfaces de fluxo vertical - North/South): coordenadas (x_centros, y_faces)
        self.X_face_y, self.Y_face_y = np.meshgrid(self._x_centros, self._y_faces, indexing='ij')

        # 3. Centros das células (para avaliação geral e plotagem): coordenadas (x_centros, y_centros)
        self.X_c, self.Y_c = np.meshgrid(self._x_centros, self._y_centros, indexing='ij')

    def _interp_2d_with_fallback(self, values_src, x_query, y_query):
        """
        Interpola um campo escalar para pontos de consulta 2D com fallback para o vizinho
        mais próximo (nearest) nas bordas (onde o método linear/cubic resulta em NaN).
        """
        # Interpolação principal (Ex: linear ou cubic)
        grid_val = griddata(self._points_orig, values_src, (x_query, y_query), method=self.interp_method)

        # Identificação e correção de pontos NaN localizados fora do fecho convexo dos pontos originais
        nan_mask = np.isnan(grid_val)
        if np.any(nan_mask):
            grid_val_nearest = griddata(self._points_orig, values_src, (x_query, y_query), method='nearest')
            grid_val[nan_mask] = grid_val_nearest[nan_mask]

        return grid_val

    def _interpolate_velocities(self):
        """
        Interpola as velocidades originais para as coordenadas das faces e centros da malha.
        """
        # Interpolação de u e v nas faces normais a x (interfaces Leste/Oeste)
        # Tamanho esperado: (Nx + 1, Ny)
        self._u_faces_x = self._interp_2d_with_fallback(self._u_orig, self.X_face_x, self.Y_face_x)
        self._v_faces_x = self._interp_2d_with_fallback(self._v_orig, self.X_face_x, self.Y_face_x)

        # Interpolação de u e v nas faces normais a y (interfaces Norte/Sul)
        # Tamanho esperado: (Nx, Ny + 1)
        self._u_faces_y = self._interp_2d_with_fallback(self._u_orig, self.X_face_y, self.Y_face_y)
        self._v_faces_y = self._interp_2d_with_fallback(self._v_orig, self.X_face_y, self.Y_face_y)

        # Interpolação de u e v nos centros das células (para gráficos e termos fonte)
        # Tamanho esperado: (Nx, Ny)
        self._u_centros = self._interp_2d_with_fallback(self._u_orig, self.X_c, self.Y_c)
        self._v_centros = self._interp_2d_with_fallback(self._v_orig, self.X_c, self.Y_c)

    # Propriedades/Atributos expostos da classe
    @property
    def x_faces(self) -> np.ndarray:
        """Coordenadas 1D das faces na direção x (tamanho: Nx + 1)"""
        return self._x_faces

    @property
    def y_faces(self) -> np.ndarray:
        """Coordenadas 1D das faces na direção y (tamanho: Ny + 1)"""
        return self._y_faces

    @property
    def x_centros(self) -> np.ndarray:
        """Coordenadas 1D dos centros na direção x (tamanho: Nx)"""
        return self._x_centros

    @property
    def y_centros(self) -> np.ndarray:
        """Coordenadas 1D dos centros na direção y (tamanho: Ny)"""
        return self._y_centros

    @property
    def u_faces(self) -> np.ndarray:
        """Velocidade u normal às faces na direção x (tamanho: Nx + 1, Ny)"""
        return self._u_faces_x

    @property
    def v_faces(self) -> np.ndarray:
        """Velocidade v normal às faces na direção y (tamanho: Nx, Ny + 1)"""
        return self._v_faces_y

    @property
    def u_centros(self) -> np.ndarray:
        """Velocidade u nos centros das células (tamanho: Nx, Ny)"""
        return self._u_centros

    @property
    def v_centros(self) -> np.ndarray:
        """Velocidade v nos centros das células (tamanho: Nx, Ny)"""
        return self._v_centros

    def plotar_vetores_velocidade(self, step: int = 1, save_path: str = None):
        """
        Visualiza o campo de velocidades interpolado nos centros das células
        usando matplotlib.pyplot.quiver e desenha as linhas que representam a malha de volumes finitos.

        Parameters
        ----------
        step : int, opcional
            Passo de amostragem dos vetores para melhorar a visibilidade do quiver. O padrão é 1.
        save_path : str, opcional
            Se fornecido, salva o gráfico gerado no caminho especificado.
        """
        fig, ax = plt.subplots(figsize=(8, 8))

        # Desenhar as linhas das faces da malha
        for x in self._x_faces:
            ax.axvline(x, color='#d3d3d3', linestyle='-', linewidth=0.5)
        for y in self._y_faces:
            ax.axhline(y, color='#d3d3d3', linestyle='-', linewidth=0.5)

        # Plotagem dos vetores
        # fatiando com step se necessário
        X = self.X_c[::step, ::step]
        Y = self.Y_c[::step, ::step]
        U = self._u_centros[::step, ::step]
        V = self._v_centros[::step, ::step]

        # Quiver
        q = ax.quiver(X, Y, U, V, color='blue', scale=5.0, width=0.003, headwidth=4)
        ax.quiverkey(q, X=0.9, Y=1.03, U=0.5, label='0.5 m/s', labelpos='E')

        # Configurações do gráfico
        ax.set_title(f"Campo de Velocidades Interpolado (Malha {self.Nx}x{self.Ny})", fontsize=14, fontweight='bold')
        ax.set_xlabel("Coordenada X (m)", fontsize=12)
        ax.set_ylabel("Coordenada Y (m)", fontsize=12)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
        ax.grid(False)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Gráfico do campo de velocidades salvo em: {save_path}")
        
        plt.show()
