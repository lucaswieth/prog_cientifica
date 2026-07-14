# Relatório Técnico de Validação de CFD
**Análise de Transporte de Escalares Acoplados por Reação (Re = 1000)**

## 1. Sumário Executivo
A auditoria presente analisa os resultados com a configuração de alta resolução espacial (**malha 150x150**) e o restabelecimento da física original (**$\Gamma = 0.001$**). O comportamento macroscópico e termodinâmico dos escalares continua fisicamente válido e o solver temporal (Euler Implícito) operou com total estabilidade. O ganho de resolução reduziu dramaticamente a "falsa difusão" do esquema UDS. No entanto, por utilizarmos $\Gamma = 0.001$, o refinamento de malha alcançado ainda não foi suficiente para trazer o **Péclet celular máximo para menos de 2**. Como o limite de estabilidade é rompido ($Pe \approx 6.67$), o esquema CDS (Diferenças Centrais) ainda induz oscilações espúrias (wiggles) nas zonas de altos gradientes convectivos, embora essas oscilações sejam muito menores em amplitude do que na malha original 40x40.

## 2. Contextualização Física e Numérica

- **Consistência Física ($\Gamma = 0.001$):** 
  O fluido simulado possui propriedades altamente convectivas (pouco difusivas). O regresso de $\Gamma$ para 0.001 exige que a malha reproduza frentes de concentração afiadas/abruptas. Conforme a taxa de reação $k$ sobe (0.0 $\rightarrow$ 0.5 $\rightarrow$ 2.0), o decaimento de $\phi_2$ é acentuado, consumindo-o em curtas distâncias de penetração no vórtice, enquanto $\phi_1$ preenche o domínio. Essa física é agora mapeada em alta definição.

- **A Discussão do Péclet (Por que o warning voltou?):**
  Ao retornar para o fluido original $\Gamma = 0.001$, a força da convecção impôs seu desafio novamente.
  Para a malha $150 \times 150$, temos $\Delta x \approx 0.00667$.
  O Número de Péclet celular máximo na tampa (onde a velocidade atinge $U \approx 1$) é matematicamente:
  $$Pe_{max} = \frac{1.0 \times 0.00667}{0.001} \approx 6.67$$
  O terminal emite o aviso *"Local Peclet number (2.52) exceeds critical value"* porque o código varre o domínio e acusa o primeiro valor que cruza a linha de corte ($Pe > 2.0$). Como $Pe_{max} \approx 6.67$ é substancialmente maior que 2, o **Critério de Scarborough continua sendo violado**.
  - **CDS (Diferenças Centrais):** Ao violar a dominância diagonal, as equações algébricas produzem coeficientes negativos. Embora o refinamento (Pe caindo de 25 para 6.67) tenha encurtado o "comprimento" das oscilações espaciais, os "wiggles" seguem vivos. O CDS ainda flutuará falsamente em regiões de bordo/camada limite.
  - **UDS (Upwind):** Este é o grande vencedor da atual configuração. Estável em qualquer Péclet, seu único defeito era a "falsa difusão", que é diretamente proporcional ao tamanho da célula ($\Delta x$). Como encolhemos $\Delta x$ para apenas $0.00667$, o UDS agora gera simulações extremamente nítidas que refletem razoavelmente bem a física de um fluido pouco difusivo (falsa difusão $\approx \Gamma$).

## 3. Análise Técnica dos Dados

- **Robustez do Solver Temporal:** 
  Avançar para uma malha de $150 \times 150 \times 2$ escalares significa montar e solucionar enormes matrizes locais a cada passo do tempo até $t=8.0s$. O fato do sistema finalizar sem nenhum erro do `scipy` ou estouro de limite (NaN) atesta perfeitamente o poder e estabilidade do **Método de Euler Implícito**. Se fosse Explícito, o limite restritivo de Courant-Friedrichs-Lewy (CFL) teria colapsado o solver com $dt = 0.1$.
  
- **Altíssima Resolução da Camada Limite:** 
  A malha agora provê pelo menos de 4 a 5 células puramente alojadas dentro da camada limite hidrodinâmica (estimada em $1/\sqrt{Re} = 0.031$). Isso traz as dinâmicas termodinâmicas capturadas ao status de altíssima fidelidade.

## 4. Validação e Conclusão Final

**Os Resultados São Válidos?**
- Para o **UDS**, sim. O comportamento atualizado une estabilidade total a uma contaminação aceitavelmente baixa por falsa difusão. É o perfil mais seguro a ser considerado real para fins práticos.
- Para o **CDS**, parcialmente. Os cortes e campos de cores continuarão demonstrando estrias não-físicas onde a velocidade forçada e o cruzamento oblíquo rasgam a diagonal matemática.

### Recomendações Definitivas
A configuração 150x150 provou que é inviável, em termos práticos, mitigar as instabilidades do CDS por puro "forçamento bruto" de malha (precisaríamos de uma malha 500x500 para ter Péclet $< 2$). A lição computacional fica evidente:

1. **Abandono do CDS Estrito:** CDS nunca deve ser a escolha principal para convecção forte, servindo apenas como bloco constitutivo acadêmico.
2. **TVD (Total Variation Diminishing):** A solução definitiva e mandatória para escoamentos $Pe \gg 2$ é incorporar um interpolador Upwind de segunda-ordem acoplado a Limitadores de Fluxo (SMART, MUSCL, QUICK), capturando frentes afiadas sem ruidos. 
3. **Malha Não-Uniforme:** Caso o uso do CDS seja imposto por restrição técnica, a expansão de células baseada em tanh ou seno garantirá densidade (apenas nas paredes) que derrubaria o Pe local, evitando o custo absurdo de uma malha 500x500 uniforme.
