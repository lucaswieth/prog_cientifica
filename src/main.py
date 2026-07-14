import os
import sys

# Ensure src is in the python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from config import SimulationConfig
from utils import setup_logger, clean_results_dir
from data_loader import VelocityDataLoader
from malha import Malha
from boundary_conditions import BoundaryConditions
from advection_schemes import get_advection_scheme
from solver import TransportSolver
from post_processing import PostProcessor

def main():
    config = SimulationConfig()
    logger = setup_logger("CFDMain")
    
    logger.info("=== Initializing Transient CFD Simulation ===")
    
    # Setup directories
    clean_results_dir(config.results_dir, logger)
    
    # 1. Load Data
    logger.info(f"Loading velocity data from {config.data_path}")
    data_loader = VelocityDataLoader(config.data_path)
    x_orig, y_orig, u_orig, v_orig = data_loader.load_data()
    
    # 2. Build Grid
    logger.info(f"Building FVM grid with resolution {config.Nx}x{config.Ny}")
    grid = Malha(x_orig, y_orig, u_orig, v_orig, config.Nx, config.Ny)
    
    # 3. Setup Boundary Conditions
    bc_phi1 = BoundaryConditions.from_dict(config.bc_phi1)
    bc_phi2 = BoundaryConditions.from_dict(config.bc_phi2)
    
    # 4. Parametric Runs
    for k in config.k_values:
        resultados_finais_phi1 = {}
        resultados_finais_phi2 = {}
        
        for scheme_name in config.schemes:
            logger.info(f"\n--- Solving Scheme: {scheme_name} | k = {k} ---")
            
            # Instantiate strategy
            advection_scheme = get_advection_scheme(scheme_name)
            
            # Instantiate solver
            solver = TransportSolver(grid, advection_scheme, config, logger)
            
            # Solve for phi_2 (Reactant)
            logger.info("Solving for phi_2 (Reactant)...")
            historico_phi2 = solver.solve_transient(
                bc_phi2, t_final=config.t_final, dt=config.dt, k=k, scalar_type='phi2'
            )
            resultados_finais_phi2[scheme_name] = historico_phi2[-1]
            
            # Post-processing phi_2
            phi2_title = f"Escalar $\\phi_2$ ({scheme_name}, k={k})"
            PostProcessor.plotar_resultados(
                historico_phi2[-1], grid, phi2_title, 
                save_path=os.path.join(config.results_dir, f"phi2_final_{scheme_name.lower()}_k{k}.png")
            )
            logger.info("Generating animation for phi_2...")
            PostProcessor.gerar_animacao_temporal(
                historico_phi2, grid, config.dt, phi2_title, 
                save_path=os.path.join(config.results_dir, f"phi2_anim_{scheme_name.lower()}_k{k}.gif")
            )
            
            # Solve for phi_1 (Product)
            logger.info("Solving for phi_1 (Product)...")
            historico_phi1 = solver.solve_transient(
                bc_phi1, t_final=config.t_final, dt=config.dt, k=k, scalar_type='phi1', phi2_hist=historico_phi2
            )
            resultados_finais_phi1[scheme_name] = historico_phi1[-1]
            
            # Post-processing phi_1
            phi1_title = f"Escalar $\\phi_1$ ({scheme_name}, k={k})"
            PostProcessor.plotar_resultados(
                historico_phi1[-1], grid, phi1_title, 
                save_path=os.path.join(config.results_dir, f"phi1_final_{scheme_name.lower()}_k{k}.png")
            )
            logger.info("Generating animation for phi_1...")
            PostProcessor.gerar_animacao_temporal(
                historico_phi1, grid, config.dt, phi1_title, 
                save_path=os.path.join(config.results_dir, f"phi1_anim_{scheme_name.lower()}_k{k}.gif")
            )
            
        # 1D Cuts (UDS vs CDS)
        logger.info(f"Generating 1D comparative cuts for k = {k}...")
        PostProcessor.plotar_corte_1d(
            resultados_finais_phi2['UDS'], resultados_finais_phi2['CDS'], grid,
            title=f"Corte em X=0.5: Escalar $\\phi_2$ (k={k})",
            save_path=os.path.join(config.results_dir, f"corte_1d_phi2_k{k}.png")
        )
        PostProcessor.plotar_corte_1d(
            resultados_finais_phi1['UDS'], resultados_finais_phi1['CDS'], grid,
            title=f"Corte em X=0.5: Escalar $\\phi_1$ (k={k})",
            save_path=os.path.join(config.results_dir, f"corte_1d_phi1_k{k}.png")
        )
        
    logger.info("=== Simulations completed successfully! ===")

if __name__ == "__main__":
    main()