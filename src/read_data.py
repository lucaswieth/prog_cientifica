import pandas as pd
import numpy as np
import os

def load_cavity_data(file_path):
    """
    Carrega o arquivo LidDrivenCavityFlow_Re_1000.dat e extrai as colunas
    x, y, u, v como arrays numpy.
    """
    # Lê o arquivo CSV
    df = pd.read_csv(file_path)
    
    # Remove espaços em branco extras dos nomes das colunas
    df.columns = df.columns.str.strip()
    
    # Extrai as colunas como arrays numpy
    x = df['x-coordinate'].to_numpy()
    y = df['y-coordinate'].to_numpy()
    u = df['x-velocity'].to_numpy()
    v = df['y-velocity'].to_numpy()
    
    return x, y, u, v

if __name__ == '__main__':
    # Caminho do arquivo de dados relativo a este script
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, 'data', 'LidDrivenCavityFlow_Re_1000.dat')
    
    print(f"Lendo o arquivo: {file_path} ...")
    
    try:
        x, y, u, v = load_cavity_data(file_path)
        
        print("\n=== Verificação dos Tipos de Dados (Numpy Arrays) ===")
        print(f"Array x: shape={x.shape}, dtype={x.dtype}")
        print(f"Array y: shape={y.shape}, dtype={y.dtype}")
        print(f"Array u: shape={u.shape}, dtype={u.dtype}")
        print(f"Array v: shape={v.shape}, dtype={v.dtype}")
        
        # Verificação se os tipos são numéricos (float) e não texto (object/string)
        if np.issubdtype(x.dtype, np.number):
            print("\n[OK] Os dados foram carregados corretamente como números (float)!")
        else:
            print("\n[AVISO] Os dados vieram formatados como texto/objeto!")
            
        print("\nAmostra dos dados (primeiras 5 linhas):")
        for i in range(5):
            print(f"Ponto {i+1}: x={x[i]:.6f}, y={y[i]:.6f}, u={u[i]:.6e}, v={v[i]:.6e}")
            
    except Exception as e:
        print(f"Erro ao ler ou processar o arquivo: {e}")
