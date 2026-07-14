# Lid-Driven Cavity Flow: Scalar Transport & Reaction

Este repositório contém uma simulação computacional de fluidos (CFD) bidimensional transiente desenvolvida em Python. O código resolve as equações de Advecção-Difusão-Reação para o transporte de dois escalares acoplados ($\phi_1$ e $\phi_2$) em uma cavidade com tampa deslizante (*Lid-Driven Cavity*), assumindo um campo de velocidades congelado (*frozen flow*) para $Re = 1000$.

O projeto utiliza o Método dos Volumes Finitos (FVM) para discretização espacial e o método de Euler Implícito para integração temporal.

## Descrição Física do Problema

A simulação analisa o transporte de dois escalares submetidos a advecção, difusão e uma reação de consumo/produção:

$$ \frac{\partial \phi_1}{\partial t} + \nabla \cdot (\vec{u} \phi_1) = \nabla \cdot (\Gamma \nabla \phi_1) + k \phi_2 $$
$$ \frac{\partial \phi_2}{\partial t} + \nabla \cdot (\vec{u} \phi_2) = \nabla \cdot (\Gamma \nabla \phi_2) - k \phi_2 $$

Onde:
- $\phi_2$ é um reagente que entra pela parede superior (Norte).
- $\phi_1$ é o produto gerado pela reação termoquímica proporcional ao consumo de $\phi_2$.
- $\Gamma$ é o coeficiente de difusão constante ($\Gamma = 0.001$).
- $k$ é a taxa de reação.
- $\vec{u}$ é o campo de velocidades (fornecido através de um dataset).

### Esquemas Numéricos

Foram implementados dois esquemas advectivos (via Padrão *Strategy*):
- **UDS (Upwind Differencing Scheme):** Incondicionalmente estável, mas de 1ª ordem e com alta difusão numérica.
- **CDS (Central Differencing Scheme):** De 2ª ordem de precisão espacial, mas suscetível a oscilações não-físicas se o critério de Péclet local superar o limite de estabilidade ($Pe = \frac{u \Delta x}{\Gamma} \le 2$). O solver está configurado para emitir avisos automatizados (warnings) se esse critério for violado.

## Estrutura do Projeto (Clean Architecture)

- `src/config.py`: Parâmetros físicos e da malha (substitui magic numbers soltos).
- `src/main.py`: Ponto de entrada (orquestrador).
- `src/solver.py`: Resolvedor das equações FVM e montagem da matriz de coeficientes lineares.
- `src/malha.py`: Gerador da geometria da malha de volumes finitos (Cell-Centered) e interpolador.
- `src/advection_schemes.py`: Implementação dos esquemas UDS e CDS.
- `src/boundary_conditions.py`: Gestão de condições de contorno de Dirichlet.
- `src/data_loader.py`: Carregamento e limpeza do campo de escoamento base via Pandas.
- `src/post_processing.py`: Visualizações e geração de animações e gráficos.
- `src/utils.py`: Funções para limpeza dos diretórios e logs.

## Instalação

Recomenda-se o uso do Anaconda/Miniconda para facilitar o gerenciamento de dependências.

```bash
# 1. Clone o repositório
git clone https://github.com/lucaswieth/prog_cientifica.git
cd prog_cientifica

# 2. Crie o ambiente virtual através do arquivo fornecido
conda env create -f environment.yml

# 3. Ative o ambiente
conda activate prog_cientifica
```

## Como Executar

Para iniciar a bateria completa de simulações paramétricas ($k = 0.0, 0.5, 2.0$ para `UDS` e `CDS`), basta executar:

```bash
python src/main.py
```

### Resultados

A pasta `results/` será limpa automaticamente em cada nova execução para evitar sobreposições. Nela, o código salvará:
- **Gráficos 2D estáticos (.png)** dos campos finais de $\phi_1$ e $\phi_2$.
- **Animações temporais (.gif)** mostrando a evolução espacial e temporal.
- **Cortes Lineares 1D (.png)** permitindo comparar a difusão numérica UDS vs CDS.

## Observações Adicionais

Todos os métodos foram refatorados respeitando conceitos de **Clean Code** e encapsulamento em **Programação Orientada a Objetos (POO)** para garantir fácil manutenção, baixa complexidade ciclomática, e escalabilidade profissional ao código de engenharia térmica.
