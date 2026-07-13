import os
import sys

# Garante que a pasta src esteja no sys.path para importações locais limpas
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from malha import Malha
from solver import TransportSolver

def main():
    # Caminhos de arquivos
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(base_dir, 'data', 'LidDrivenCavityFlow_Re_1000.dat')
    artifact_dir = os.path.join(base_dir, 'results')
    os.makedirs(artifact_dir, exist_ok=True)

    print("=== Inicializando Simulação CFD Transiente ===")
    Nx, Ny = 40, 40
    print(f"Resolução da malha: {Nx}x{Ny} volumes finitos")
    
    # 1. Carrega e gera a malha
    grid = Malha(data_path, Nx=Nx, Ny=Ny, interp_method='linear')
    
    # 2. Configura as condições de contorno (Dirichlet)
    bc_phi1 = {
        'west': 0.0,
        'east': 0.0,
        'south': 1.0,
        'north': 0.0
    }
    
    bc_phi2 = {
        'west': 0.0,
        'east': 0.0,
        'south': 0.0,
        'north': 1.0
    }
    
    # 3. Bateria de Simulações Paramétricas
    schemes = ['UDS', 'CDS']
    k_values = [0.0, 0.5, 2.0]
    
    # Parâmetros Temporais
    t_final = 8.0
    dt = 0.1
    
    for k in k_values:
        # Armazenaremos os resultados finais para plotar o corte 1D comparativo
        resultados_finais_phi1 = {}
        resultados_finais_phi2 = {}
        
        for scheme in schemes:
            print(f"\n--- Resolvendo com o esquema: {scheme} | k = {k} ---")
            
            # Instancia o solucionador
            solver = TransportSolver(grid, scheme=scheme, gamma=0.001)
            
            # 1. Resolve para phi_2 (Reagente consumido, fonte na parede superior)
            print("  Resolvendo para phi_2 (Reagente) no tempo...")
            historico_phi2 = solver.resolver_transiente(bc_phi2, t_final=t_final, dt=dt, k=k, scalar_type='phi2')
            resultados_finais_phi2[scheme] = historico_phi2[-1]
            
            # Gráficos finais e Animação phi_2
            phi2_title = f"Escalar $\\phi_2$ ({scheme}, k={k})"
            solver.plotar_resultados(historico_phi2[-1], phi2_title, save_path=os.path.join(artifact_dir, f"phi2_final_{scheme.lower()}_k{k}.png"))
            print("  Gerando animação para phi_2...")
            solver.gerar_animacao_temporal(historico_phi2, dt, phi2_title, save_path=os.path.join(artifact_dir, f"phi2_anim_{scheme.lower()}_k{k}.gif"))
            
            # 2. Resolve para phi_1 (Produto gerado, fonte na parede inferior)
            print("  Resolvendo para phi_1 (Produto, usando phi_2 como fonte)...")
            historico_phi1 = solver.resolver_transiente(bc_phi1, t_final=t_final, dt=dt, k=k, scalar_type='phi1', phi2_hist=historico_phi2)
            resultados_finais_phi1[scheme] = historico_phi1[-1]
            
            # Gráficos finais e Animação phi_1
            phi1_title = f"Escalar $\\phi_1$ ({scheme}, k={k})"
            solver.plotar_resultados(historico_phi1[-1], phi1_title, save_path=os.path.join(artifact_dir, f"phi1_final_{scheme.lower()}_k{k}.png"))
            print("  Gerando animação para phi_1...")
            solver.gerar_animacao_temporal(historico_phi1, dt, phi1_title, save_path=os.path.join(artifact_dir, f"phi1_anim_{scheme.lower()}_k{k}.gif"))
            
        # 3. Cortes 1D Comparativos (UDS vs CDS) para o k atual
        print(f"\n  Gerando cortes 1D comparativos para k = {k}...")
        TransportSolver.plotar_corte_1d(
            phi_uds=resultados_finais_phi2['UDS'], 
            phi_cds=resultados_finais_phi2['CDS'], 
            grid=grid, 
            title=f"Corte em X=0.5: Escalar $\\phi_2$ (k={k})", 
            save_path=os.path.join(artifact_dir, f"corte_1d_phi2_k{k}.png")
        )
        TransportSolver.plotar_corte_1d(
            phi_uds=resultados_finais_phi1['UDS'], 
            phi_cds=resultados_finais_phi1['CDS'], 
            grid=grid, 
            title=f"Corte em X=0.5: Escalar $\\phi_1$ (k={k})", 
            save_path=os.path.join(artifact_dir, f"corte_1d_phi1_k{k}.png")
        )

    print("\n=== Simulações concluídas com sucesso! ===")

if __name__ == '__main__':
    main()