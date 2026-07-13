import os
import sys
import numpy as np
import pytest

# Adiciona o diretório src ao sys.path para permitir a importação dos módulos do projeto
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from malha import Malha
from solver import TransportSolver

@pytest.fixture
def grid_inst():
    """Fixture para inicializar uma instância da classe Malha de tamanho 40x40."""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(base_dir, 'data', 'LidDrivenCavityFlow_Re_1000.dat')
    return Malha(data_path, Nx=40, Ny=40, interp_method='linear')

def test_solver_creation(grid_inst):
    """Verifica se o solver é criado corretamente e rejeita esquemas inválidos."""
    solver = TransportSolver(grid_inst, scheme='UDS', gamma=1.0)
    assert solver.scheme == 'UDS'
    assert solver.gamma == 1.0
    
    with pytest.raises(ValueError):
        TransportSolver(grid_inst, scheme='INVALID_SCHEME')

def test_solver_resolution_uds(grid_inst):
    """Valida o resultado do solver com o esquema UDS (deve respeitar os limites de Dirichlet)."""
    solver = TransportSolver(grid_inst, scheme='UDS', gamma=1.0)
    
    # Condições de contorno para phi_1 (fonte na parede inferior)
    bc_phi1 = {
        'west': 0.0,
        'east': 0.0,
        'south': 1.0,
        'north': 0.0
    }
    
    phi1_2d = solver.resolver(bc_phi1)
    
    # Verifica o formato e se não existem valores nulos (NaN)
    assert phi1_2d.shape == (40, 40)
    assert not np.isnan(phi1_2d).any()
    
    # Sob o esquema UDS, a solução deve respeitar estritamente o princípio do máximo:
    # os valores devem estar limitados entre os extremos das condições de contorno [0.0, 1.0].
    assert phi1_2d.min() >= -1e-12
    assert phi1_2d.max() <= 1.0 + 1e-12

def test_solver_resolution_cds(grid_inst):
    """Valida o resultado do solver com o esquema CDS."""
    solver = TransportSolver(grid_inst, scheme='CDS', gamma=1.0)
    
    # Condições de contorno para phi_2 (fonte na parede superior)
    bc_phi2 = {
        'west': 0.0,
        'east': 0.0,
        'south': 0.0,
        'north': 1.0
    }
    
    phi2_2d = solver.resolver(bc_phi2)
    
    # Verifica o formato e se não existem valores nulos (NaN)
    assert phi2_2d.shape == (40, 40)
    assert not np.isnan(phi2_2d).any()
