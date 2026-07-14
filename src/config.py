from dataclasses import dataclass, field
from typing import List
import os

@dataclass
class SimulationConfig:
    """
    Centralizes all physical and numerical parameters for the CFD simulation.
    """
    
    # Mesh resolution
    Nx: int = 150
    Ny: int = 150
    
    # Physical parameters
    gamma: float = 0.001
    k_values: List[float] = field(default_factory=lambda: [0.0, 0.5, 2.0])
    
    # Time integration parameters
    t_final: float = 8.0
    dt: float = 0.1
    
    # Numerical schemes to evaluate
    schemes: List[str] = field(default_factory=lambda: ['UDS', 'CDS'])
    
    # Stability criterion
    critical_peclet: float = 2.0
    
    # Boundary conditions definitions
    bc_phi1: dict = field(default_factory=lambda: {
        'west': 0.0,
        'east': 0.0,
        'south': 1.0,
        'north': 0.0
    })
    
    bc_phi2: dict = field(default_factory=lambda: {
        'west': 0.0,
        'east': 0.0,
        'south': 0.0,
        'north': 1.0
    })

    # File paths
    base_dir: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path: str = os.path.join(base_dir, 'data', 'LidDrivenCavityFlow_Re_1000.dat')
    results_dir: str = os.path.join(base_dir, 'results')
