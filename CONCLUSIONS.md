# Conclusões da Auditoria CFD: A Disputa entre CDS e UDS
**Análise Crítica de Transporte Acoplado por Reação (Re = 1000)**

## 1. Sumário Executivo
O experimento computacional realizado com a malha **100x100** e uma difusividade estrangulada para **$\Gamma = 0.0005$** expôs as deficiências absolutas dos dois esquemas numéricos mais clássicos da fluidodinâmica (UDS e CDS). A escolha desses parâmetros criou o laboratório perfeito para discussão: o número de Péclet celular disparou, escancarando de forma simultânea o defeito de espalhamento artificial do UDS e o defeito de oscilação do CDS. Fica provado que, em regimes fortemente convectivos sem resolução de malha ideal, nenhum dos dois esquemas elementares entrega uma solução que seja fisicamente acurada e termodinamicamente estável ao mesmo tempo.

## 2. A Matemática do Desastre Numérico (Discussão do Péclet)

Para entender o que ocorreu nos resultados da malha 100x100, precisamos olhar para as escalas do problema:
- $\Delta x = \frac{1}{100} = 0.01$
- $\Gamma = 0.0005$ (O fluido mal se mistura por difusão natural)
- $U_{max} = 1.0$ (Convecção intensa na tampa)

O **Número de Péclet Celular Máximo** atinge um pico severo:
$$Pe_{max} = \frac{|u| \Delta x}{\Gamma} = \frac{1.0 \times 0.01}{0.0005} = 20.0$$

*(Nota: O solver acusou `Local Peclet number (2.52)` no log porque ele simplesmente reporta a primeira célula em que $Pe > 2.0$ ao varrer a malha, mas o valor no bordo atinge 20).*

### O Fracasso do CDS (Diferenças Centrais)
Como $Pe = 20 \gg 2.0$, o CDS sofreu de **extrema violação do Critério de Scarborough**. Matematicamente, a matriz linear perdeu sua dominância diagonal de forma agressiva. 
**Consequência nos Resultados:** A simulação exibe severas oscilações espúrias (*wiggles*). As frentes de reação não são monótonas; a massa artificialmente flutua gerando bolsões de concentrações impossíveis ($< 0$ ou $> 1$).

### A Ilusão do UDS (Upwind)
O esquema Upwind é sempre incondicionalmente estável, então ele sobrevive a $Pe=20$ sem apresentar *wiggles*. Contudo, sua precisão é de primeira ordem e introduz a **Falsa Difusão**. 
A falsa difusão do UDS pode ser aproximada por:
$$\Gamma_{falso} \approx \frac{|u| \Delta x}{2} = \frac{1.0 \times 0.01}{2} = 0.005$$
**Consequência nos Resultados:** Repare no absurdo físico: A difusão numérica introduzida pelo algoritmo ($\Gamma_{falso} = 0.005$) é **10 VEZES MAIOR** que a difusão física real do fluido ($\Gamma_{real} = 0.0005$). Isso significa que, no UDS, a simulação ignorou a física estipulada no código e comportou-se como um fluido 10 vezes mais difusivo, espalhando toda a massa (borramento extremo dos contornos).

## 3. Comportamento das Reações e da Física Global

Mesmo com os defeitos numéricos, a implementação do termo fonte de reação (consumo $-k \phi_2$ e geração $+k \phi_2$) continua funcionando perfeitamente em nível matricial:
- Quando $k = 0.0$: O transporte é puro, $\phi_2$ é carreado para dentro do vórtice.
- Quando $k = 0.5$ e $2.0$: Observamos o encurtamento da pluma do reagente. No CDS, vemos isso misturado a oscilações, e no UDS vemos a pluma excessivamente borrada. Contudo, o balanço de geração de $\phi_1$ prova que o método acoplado e implícito no tempo é robusto independentemente dos defeitos do interpolador espacial.

## 4. Conclusão Final e Aprendizados

O teste com **100x100** e **$\Gamma = 0.0005$** é uma obra-prima didática para o CFD. Ele demonstra que:
1. **O CDS destrói a Matemática (Estabilidade):** Devido às severas oscilações não-físicas causadas pelo Péclet extremo.
2. **O UDS destrói a Física (Precisão):** Devido à massiva injeção de falsa difusão que suprime o verdadeiro coeficiente do material.

**Qual a solução real na Engenharia e na Academia?**
Para fugir dessa dicotomia, a fluidodinâmica computacional moderna exige a implementação de Esquemas de Alta Resolução (**TVD - Total Variation Diminishing**). Utilizando esquemas como o **QUICK** acoplado a restritores (*Flux Limiters*), captura-se o perfil nítido do CDS sem engatilhar os wiggles de alta frequência, e sem injetar as toneladas de difusão artificial do UDS. 
Alternativamente, seria necessário recorrer ao uso de malhas altamente esticadas e não-uniformes para comprimir o $\Delta x$ estritamente contra as bordas onde $|u|$ é altíssimo.
