import os
import sys
import numpy as np
import pytest

# Adiciona o diretório src ao sys.path para permitir a importação dos módulos do projeto
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from malha import Malha

@pytest.fixture
def grid_inst():
    """Fixture para inicializar uma instância da classe Malha de tamanho 40x40."""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(base_dir, 'data', 'LidDrivenCavityFlow_Re_1000.dat')
    return Malha(data_path, Nx=40, Ny=40, interp_method='linear')

def test_geometria_malha(grid_inst):
    """Verifica as dimensões físicas e de discretização da malha."""
    Nx, Ny = 40, 40
    assert len(grid_inst.x_faces) == Nx + 1
    assert len(grid_inst.y_faces) == Ny + 1
    assert len(grid_inst.x_centros) == Nx
    assert len(grid_inst.y_centros) == Ny
    
    # Domínio físico deve ser [0, 1] x [0, 1]
    assert np.isclose(grid_inst.x_faces.min(), 0.0)
    assert np.isclose(grid_inst.x_faces.max(), 1.0)
    assert np.isclose(grid_inst.y_faces.min(), 0.0)
    assert np.isclose(grid_inst.y_faces.max(), 1.0)

def test_valores_centros(grid_inst):
    """Garante que os centros sejam as médias das faces vizinhas."""
    assert np.isclose(grid_inst.x_centros[0], 0.5 * (grid_inst.x_faces[0] + grid_inst.x_faces[1]))
    assert np.isclose(grid_inst.y_centros[0], 0.5 * (grid_inst.y_faces[0] + grid_inst.y_faces[1]))

def test_velocidades_malha(grid_inst):
    """Valida o formato e integridade das velocidades interpoladas."""
    Nx, Ny = 40, 40
    assert grid_inst.u_faces.shape == (Nx + 1, Ny)
    assert grid_inst.v_faces.shape == (Nx, Ny + 1)
    assert grid_inst.u_centros.shape == (Nx, Ny)
    assert grid_inst.v_centros.shape == (Nx, Ny)
    
    # Garante que não existem valores nulos (NaN)
    assert not np.isnan(grid_inst.u_faces).any(), "Existem NaNs em u_faces!"
    assert not np.isnan(grid_inst.v_faces).any(), "Existem NaNs em v_faces!"
    assert not np.isnan(grid_inst.u_centros).any(), "Existem NaNs em u_centros!"
    assert not np.isnan(grid_inst.v_centros).any(), "Existem NaNs em v_centros!"
