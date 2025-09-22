# Modelo de Reator Químico com Runge-Kutta de 4ª Ordem

**Disciplina:** Métodos Numéricos para Equações Diferenciais (MNED)
**Instituição:** Instituto Politécnico do Rio de Janeiro (IPRJ/UERJ)
**Autores:**
* Lucas Rodrigues Estorck Pinto (202310349511)
* Tainá Martins Macário (202310070411)
**Professor:** Grazione Souza
**Data:** 20 de setembro de 2025

---

## 1. Visão Geral do Projeto

Este projeto consiste na implementação do método de Runge-Kutta clássico de quarta ordem (RK4) para a solução numérica de um sistema de equações diferenciais ordinárias (EDOs). O sistema modela a dinâmica de um reator de batelada para a produção de penicilina por fermentação.

O trabalho foi desenvolvido na linguagem de programação Julia, escolhida por seu alto desempenho em computação científica.

## 2. O Modelo Matemático

O comportamento do reator é descrito pelo seguinte sistema de EDOs:

$$
\frac{dg}{dt} = 13.1g - 13.94g^2
$$

$$
\frac{df}{dt} = 1.71g
$$

Onde:
* `g`: Concentração adimensional da massa das células.
* `f`: Concentração adimensional de penicilina.
* `t`: Tempo adimensional no intervalo $0 \le t \le 1$.

As condições iniciais são:
* $g(0) = 0.03$
* $f(0) = 0.0$

## 3. Estrutura do Repositório

```
/
├── .gitignore          # Arquivos e pastas a serem ignorados pelo Git
├── LICENSE             # Licença MIT do projeto
├── README.md           # Documentação principal
├── assignment/         # Contém o enunciado original do trabalho
├── report/             # Contém o relatório final em PDF
└── src/                # Contém o código-fonte em Julia
```

## 4. Como Executar o Código

### Pré-requisitos
* **Julia:** É necessário ter o Julia (versão 1.11.7+ ou superior) instalado. Pode baixá-lo em [julialang.org](https://julialang.org/).

### Instalação de Dependências

O código utiliza as seguintes bibliotecas: `Plots`, `Printf`, `Statistics`, e `Dates`. Para instalá-las, abra o terminal do Julia (REPL) e execute:

```julia
using Pkg
Pkg.add(["Plots", "Printf", "Statistics", "Dates"])
```

### Execução

1.  Navegue até à pasta raiz do projeto no seu terminal.
2.  Execute o script principal com o comando:

```bash
julia src/reactor_simulation.jl
```

Ao ser executado, o script realizará a simulação, a análise de convergência e a análise de sensibilidade, gerando os seguintes arquivos na pasta raiz:

* **Tabelas de dados (.txt):** `resultados_dt_*.txt`, `convergencia_numerica.txt`, `analise_sensibilidade_completa.txt`.
* **Gráficos (.png):** `evolucao_g_tempo.png`, `evolucao_f_tempo.png`, `plano_fase.png`, e diversos gráficos de sensibilidade.

## 5. Principais Resultados

A simulação demonstra um crescimento logístico para a concentração de células (`g(t)`) e um crescimento quase linear para a concentração de penicilina (`f(t)`) após uma fase inicial.

Para uma análise completa e discussão detalhada, consulte o [relatório técnico](report/MEDO___Relatório_1.pdf) neste repositório.

## 6. Licença

Este projeto está licenciado sob a Licença MIT. Veja o arquivo [LICENSE](LICENSE) para mais detalhes.

[EN]

# 4th-Order Runge-Kutta for a Chemical Reactor Model

**Course:** Numerical Methods for Differential Equations (MNED)
**Institution:** Polytechnic Institute of Rio de Janeiro (IPRJ/UERJ)
**Authors:**
* Lucas Rodrigues Estorck Pinto (202310349511)
* Tainá Martins Macário (202310070411)
**Professor:** Grazione Souza
**Date:** September 20, 2025

---

## 1. Project Overview

This project consists of the implementation of the classic fourth-order Runge-Kutta (RK4) method for the numerical solution of a system of ordinary differential equations (ODEs). The system models the dynamics of a batch reactor for penicillin production through fermentation.

The project was developed in the Julia programming language, chosen for its high performance in scientific computing.

## 2. The Mathematical Model

[cite_start]The reactor's behavior is described by the following system of ODEs[cite: 24, 25]:

$$
\frac{dg}{dt} = 13.1g - 13.94g^2 
$$

$$
\frac{df}{dt} = 1.71g 
$$

[cite_start]Where[cite: 27]:
* `g`: Dimensionless cell mass concentration.
* `f`: Dimensionless penicillin concentration.
* `t`: Dimensionless time in the interval $0 \le t \le 1$.

[cite_start]The initial conditions are[cite: 28, 29]:
* $g(0) = 0.03$
* $f(0) = 0.0$

## 3. Repository Structure

```
/
├── .gitignore          # Files and folders to be ignored by Git
├── LICENSE             # MIT License for the project
├── README.md           # Main documentation file (English)
├── README-pt.md        # Documentation in Portuguese
├── assignment/         # Contains the original assignment description
├── report/             # Contains the final report in PDF
└── src/                # Contains the Julia source code
```

## 4. How to Run the Code

### Prerequisites
* **Julia:** Julia (version 1.11.7+ or higher) must be installed. You can download it from [julialang.org](https://julialang.org/).

### Installing Dependencies

The code uses the following libraries: `Plots`, `Printf`, `Statistics`, and `Dates`. To install them, open the Julia REPL (terminal) and run:

```julia
using Pkg
Pkg.add(["Plots", "Printf", "Statistics", "Dates"])
```

### Execution

1.  Navigate to the project's root directory in your terminal.
2.  Run the main script using the command:

```bash
julia src/reactor_simulation.jl
```

When executed, the script will perform the simulation, convergence analysis, and sensitivity analysis, generating the following files in the root directory:

* **Data tables (.txt):** `resultados_dt_*.txt`, `convergencia_numerica.txt`, `analise_sensibilidade_completa.txt`.
* **Plots (.png):** `evolucao_g_tempo.png`, `evolucao_f_tempo.png`, `plano_fase.png`, and several sensitivity analysis plots.

## 5. Key Results

[cite_start]The simulation shows a logistic growth for the cell concentration (`g(t)`) and a nearly linear growth for the penicillin concentration (`f(t)`) after an initial transient phase[cite: 345].

For a complete analysis and detailed discussion, please refer to the [technical report](report/MEDO___Relatório_1.pdf) in this repository. (Note: The report is in Portuguese).

## 6. License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
