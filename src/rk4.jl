using Plots, Printf, Statistics, Dates

# Estrutura para armazenar parâmetros do modelo
struct ReactorParams
    a1::Float64  # Coeficiente linear (13.1)
    a2::Float64  # Coeficiente quadrático (13.94)
    b1::Float64  # Coeficiente produção penicilina (1.71)
end

# Parâmetros originais conforme especificação do problema
PARAMS_ORIGINAL = ReactorParams(13.1, 13.94, 1.71)

function dg_dt(g, f, t, params::ReactorParams)
    return params.a1 * g - params.a2 * g^2
end

function df_dt(g, f, t, params::ReactorParams)
    return params.b1 * g
end

# Implementação do método Runge-Kutta clássico de 4 ordem
function runge_kutta_4(g0, f0, dt, t_max, params::ReactorParams)
    """
    Implementa o método Runge-Kutta clássico de 4ª ordem para sistema de EDOs
    
    Parâmetros:
    - g0, f0: condições iniciais
    - dt: incremento de tempo
    - t_max: tempo máximo de simulação
    - params: parâmetros do modelo
    
    Retorna:
    - t, g, f: vetores com soluções numéricas
    """
    
    # Número de passos
    n_steps = Int(ceil(t_max / dt))
    
    # Inicializar vetores de resultados
    t = Vector{Float64}(undef, n_steps + 1)
    g = Vector{Float64}(undef, n_steps + 1)
    f = Vector{Float64}(undef, n_steps + 1)
    
    # Condições iniciais
    t[1] = 0.0
    g[1] = g0
    f[1] = f0
    
    # Loop principal do método RK4
    for i in 1:n_steps
        # Valores atuais
        t_n = t[i]
        g_n = g[i]
        f_n = f[i]
        
        # Cálculo dos coeficientes k para g
        k1_g = dg_dt(g_n, f_n, t_n, params)
        k2_g = dg_dt(g_n + k1_g * dt/2, f_n, t_n + dt/2, params)
        k3_g = dg_dt(g_n + k2_g * dt/2, f_n, t_n + dt/2, params)
        k4_g = dg_dt(g_n + k3_g * dt, f_n, t_n + dt, params)
        
        # Cálculo dos coeficientes k para f
        k1_f = df_dt(g_n, f_n, t_n, params)
        k2_f = df_dt(g_n + k1_g * dt/2, f_n, t_n + dt/2, params)
        k3_f = df_dt(g_n + k2_g * dt/2, f_n, t_n + dt/2, params)
        k4_f = df_dt(g_n + k3_g * dt, f_n, t_n + dt, params)
        
        # Aplicação da fórmula RK4
        g[i+1] = g_n + (dt/6) * (k1_g + 2*k2_g + 2*k3_g + k4_g)
        f[i+1] = f_n + (dt/6) * (k1_f + 2*k2_f + 2*k3_f + k4_f)
        t[i+1] = t_n + dt
    end
    
    return t, g, f
end

# Estrutura para armazenar resultados
struct SimulationResult
    t::Vector{Float64}
    g::Vector{Float64}
    f::Vector{Float64}
    params::ReactorParams
    dt::Float64
    label::String
end

# Função para salvar resultados em formato tabular
function salvar_resultados_tabela(t, g, f, dt, filename)
    """
    Salva os resultados em formato de tabela para o relatório
    """
    open(filename, "w") do file
        write(file, "="^70 * "\n\n")
        
        write(file, "PARÂMETROS DA SIMULAÇÃO:\n")
        write(file, @sprintf("• Incremento de tempo (Δt): %.4f\n", dt))
        write(file, @sprintf("• Tempo máximo (t_max): %.1f\n", t[end]))
        write(file, @sprintf("• Número de passos: %d\n", length(t)-1))
        write(file, "• Condições iniciais: g(0) = 0.03, f(0) = 0.0\n")
        write(file, "• Coeficientes: a₁ = 13.1, a₂ = 13.94, b₁ = 1.71\n\n")
        
        write(file, "RESULTADOS NUMÉRICOS:\n")
        write(file, @sprintf("%-12s %-15s %-15s\n", "t", "g(t)", "f(t)"))
        write(file, "-"^45 * "\n")
        
        # Mostrar todos os pontos para dt principal, ou pontos selecionados para outros
        step = dt == 0.05 ? 1 : max(1, Int(round(length(t)/20)))
        
        for i in 1:step:length(t)
            write(file, @sprintf("%-12.4f %-15.8f %-15.8f\n", t[i], g[i], f[i]))
        end
        
        # Sempre mostrar o último ponto
        if step > 1
            write(file, @sprintf("%-12.4f %-15.8f %-15.8f\n", t[end], g[end], f[end]))
        end
        
        write(file, "\n")
        write(file, @sprintf("VALORES FINAIS (t = %.1f):\n", t[end]))
        write(file, @sprintf("• g(%.1f) = %.8f\n", t[end], g[end]))
        write(file, @sprintf("• f(%.1f) = %.8f\n", t[end], f[end]))

        # Análise específica para Δt altos
        if dt >= 0.5
            write(file, "\nANÁLISE PARA Δt ALTO:\n")
            write(file, "• Com Δt = %.1f, utilizamos apenas %d passo(s)\n", dt, length(t)-1)
            if dt >= 1.0
                write(file, "• ATENÇÃO: Δt ≥ t_max pode causar instabilidade numérica\n")
            else
                write(file, "• Δt alto pode comprometer a precisão da solução\n")
            end
        end
    end
    println("Resultados salvos em: $filename")
end

# Função para análise de convergência conforme especificação
function analise_convergencia_completa(g0, f0, t_max)
    """
    Realiza análise de convergência variando Δt conforme requisito 2
    Incremento inicial: Δt = 0.05, depois testa convergência
    """
    println("\n" * "="^60)
    println("ANÁLISE DE CONVERGÊNCIA NUMÉRICA")
    println("="^60)
    
    # Valores de Δt para teste de convergência
    dt_values = [1.0, 0.7, 0.05, 0.025, 0.01, 0.005, 0.001]
    resultados = SimulationResult[]
    
    println("Testando convergência com diferentes valores de Δt:")
    println("-"^60)
    
    # Solução de referência com Δt pequeno
    println("Calculando solução de referência (Δt = 0.0005)...")
    t_ref, g_ref, f_ref = runge_kutta_4(g0, f0, 0.0005, t_max, PARAMS_ORIGINAL)
    
    # Arquivo para resultados de convergência
    open("convergencia_numerica.txt", "w") do file
        write(file, "ANÁLISE DE CONVERGÊNCIA NUMÉRICA\n")
        write(file, "="^70 * "\n\n")
        
        write(file, @sprintf("%-8s %-12s %-12s %-12s %-12s %-10s\n", 
              "Δt", "g(1)", "f(1)", "Erro g", "Erro f", "Ordem"))
        write(file, "-"^75 * "\n")
        
        erro_anterior_g = 0.0
        
        for (idx, dt) in enumerate(dt_values)
            println(@sprintf("Simulando com Δt = %.4f...", dt))
            
            t, g, f = runge_kutta_4(g0, f0, dt, t_max, PARAMS_ORIGINAL)
            
            # Calcular erros em relação à solução de referência
            erro_g = abs(g[end] - g_ref[end])
            erro_f = abs(f[end] - f_ref[end])
            
            # Calcular ordem de convergência
            ordem = idx > 1 && erro_anterior_g > 0 ? log2(erro_anterior_g / erro_g) : 0.0
            
            write(file, @sprintf("%-8.4f %-12.8f %-12.8f %-12.2e %-12.2e %-10.2f\n", 
                  dt, g[end], f[end], erro_g, erro_f, ordem))
            
            println(@sprintf("  Δt = %.4f: g(1) = %.8f, f(1) = %.8f, Erro_g = %.2e", 
                    dt, g[end], f[end], erro_g))
            
            # Salvar resultados tabulares para cada Δt
            filename = @sprintf("resultados_dt_%.4f.txt", dt)
            salvar_resultados_tabela(t, g, f, dt, filename)
            
            # Armazenar para plotagem
            resultado = SimulationResult(t, g, f, PARAMS_ORIGINAL, dt, "Δt = $dt")
            push!(resultados, resultado)
            
            erro_anterior_g = erro_g
        end
        
        write(file, "\nNOTAS:\n")
        write(file, "• Solução de referência calculada com Δt = 0.0005\n")
        write(file, "• Ordem de convergência teórica do RK4: 4\n")
        write(file, "• Erro = |valor_calculado - valor_referência|\n")
    end
    
    return resultados
end

# Função para plotar resultados principais
function plotar_resultados_principais(resultados)
    """
    Gera os gráficos principais dos resultados conforme especificação
    SEPARA gráficos de Δt altos dos normais
    """
    println("\nGerando gráficos dos resultados...")
    
    # Configurar tema dos gráficos
    gr(size=(800, 600), dpi=300)
    
    # Separar resultados por Δt
    resultados_normais = filter(r -> r.dt <= 0.05, resultados)
    resultados_altos = filter(r -> r.dt > 0.05, resultados)
    
    # === GRÁFICOS NORMAIS (Δt ≤ 0.05) ===
    p1 = plot(title="Evolução da Concentração de Células g(t)", 
              xlabel="Tempo adimensional (t)", 
              ylabel="Concentração adimensional g(t)",
              legend=:bottomright,
              grid=true,
              linewidth=2)
    
    p2 = plot(title="Evolução da Concentração de Penicilina f(t)", 
              xlabel="Tempo adimensional (t)", 
              ylabel="Concentração adimensional f(t)",
              legend=:bottomright,
              grid=true,
              linewidth=2)
    
    p3 = plot(title="Plano de Fase: g(t) vs f(t)", 
              xlabel="Concentração de células g(t)", 
              ylabel="Concentração de penicilina f(t)",
              legend=:bottomright,
              grid=true,
              linewidth=2)
    
    # Cores para Δt normais
    cores_normais = [:red, :blue, :green, :orange, :purple]
    
    for (i, res) in enumerate(resultados_normais)
        cor = cores_normais[min(i, length(cores_normais))]
        estilo = res.dt == 0.05 ? :solid : :dash
        largura = res.dt == 0.05 ? 3 : 2
        
        plot!(p1, res.t, res.g, label=res.label, color=cor, linestyle=estilo, linewidth=largura)
        plot!(p2, res.t, res.f, label=res.label, color=cor, linestyle=estilo, linewidth=largura)
        plot!(p3, res.g, res.f, label=res.label, color=cor, linestyle=estilo, linewidth=largura)
    end
    
    # === GRÁFICOS ALTOS ===
    if !isempty(resultados_altos)
        p4 = plot(title="Evolução de g(t) - Δt ALTOS", 
                  xlabel="Tempo adimensional (t)", 
                  ylabel="Concentração adimensional g(t)",
                  legend=:bottomright,
                  grid=true,
                  linewidth=3)
        
        p5 = plot(title="Evolução de f(t) - Δt ALTOS", 
                  xlabel="Tempo adimensional (t)", 
                  ylabel="Concentração adimensional f(t)",
                  legend=:bottomright,
                  grid=true,
                  linewidth=3)
        
        p6 = plot(title="Plano de Fase - Δt ALTOS", 
                  xlabel="Concentração de células g(t)", 
                  ylabel="Concentração de penicilina f(t)",
                  legend=:bottomright,
                  grid=true,
                  linewidth=3)
        
        # Cores para Δt altos
        cores_altos = [:red, :orange, :purple]
        
        for (i, res) in enumerate(resultados_altos)
            cor = cores_altos[min(i, length(cores_altos))]
            marker = res.dt >= 1.0 ? :circle : :square
            markersize = res.dt >= 1.0 ? 8 : 6
            
            plot!(p4, res.t, res.g, label=res.label, color=cor, 
                  marker=marker, markersize=markersize, linewidth=3)
            plot!(p5, res.t, res.f, label=res.label, color=cor, 
                  marker=marker, markersize=markersize, linewidth=3)
            plot!(p6, res.g, res.f, label=res.label, color=cor, 
                  marker=marker, markersize=markersize, linewidth=3)
        end
        
        # Salvar gráficos de Δt altos
        savefig(p4, "evolucao_g_dt_altos.png")
        savefig(p5, "evolucao_f_dt_altos.png")
        savefig(p6, "plano_fase_dt_altos.png")
        
        # Gráfico combinado para Δt altos
        p_combined_altos = plot(p4, p5, layout=(2,1), size=(800, 1000))
        savefig(p_combined_altos, "resultados_dt_altos.png")
    end
    
    # Salvar gráficos normais
    savefig(p1, "evolucao_g_tempo.png")
    savefig(p2, "evolucao_f_tempo.png")
    savefig(p3, "plano_fase.png")
    
    # Gráfico combinado normal
    p_combined = plot(p1, p2, layout=(2,1), size=(800, 1000))
    savefig(p_combined, "resultados_principais.png")
    
    println("Gráficos principais salvos:")
    println("• evolucao_g_tempo.png (Δt normais)")
    println("• evolucao_f_tempo.png (Δt normais)")
    println("• plano_fase.png (Δt normais)")
    println("• resultados_principais.png (Δt normais)")
    
    if !isempty(resultados_altos)
        println("\nGráficos separados para Δt altos:")
        println("• evolucao_g_dt_altos.png")
        println("• evolucao_f_dt_altos.png") 
        println("• plano_fase_dt_altos.png")
        println("• resultados_dt_altos.png")
    end
    
    return p1, p2, p3
end

# Função para análise de sensibilidade
function analise_sensibilidade_coeficientes(g0, f0, t_max, dt_base)
    """
    Análise de sensibilidade quando da variação de coeficientes das equações
    Conforme requisito 4 da especificação
    """
    println("\n" * "="^60)
    println("ANÁLISE DE SENSIBILIDADE DOS COEFICIENTES")
    println("="^60)
    
    # Variações percentuais dos coeficientes
    variacoes = [-10, -5, 0, 5, 10]  # em %
    resultados_sensibilidade = Dict{String, Vector{SimulationResult}}()
    
    # Coeficientes originais
    a1_orig, a2_orig, b1_orig = 13.1, 13.94, 1.71
    
    # Análise do coeficiente a₁ = 13.1
    println("Analisando sensibilidade do coeficiente a₁ = 13.1...")
    resultados_a1 = SimulationResult[]
    
    for var in variacoes
        fator = 1 + var/100
        a1_novo = a1_orig * fator
        params_mod = ReactorParams(a1_novo, a2_orig, b1_orig)
        
        t, g, f = runge_kutta_4(g0, f0, dt_base, t_max, params_mod)
        label = var == 0 ? "Original (a₁=13.1)" : @sprintf("a₁=%.2f (%+d%%)", a1_novo, var)
        resultado = SimulationResult(t, g, f, params_mod, dt_base, label)
        push!(resultados_a1, resultado)
    end
    resultados_sensibilidade["a1"] = resultados_a1
    
    # Análise do coeficiente a₂ = 13.94
    println("Analisando sensibilidade do coeficiente a₂ = 13.94...")
    resultados_a2 = SimulationResult[]
    
    for var in variacoes
        fator = 1 + var/100
        a2_novo = a2_orig * fator
        params_mod = ReactorParams(a1_orig, a2_novo, b1_orig)
        
        t, g, f = runge_kutta_4(g0, f0, dt_base, t_max, params_mod)
        label = var == 0 ? "Original (a₂=13.94)" : @sprintf("a₂=%.2f (%+d%%)", a2_novo, var)
        resultado = SimulationResult(t, g, f, params_mod, dt_base, label)
        push!(resultados_a2, resultado)
    end
    resultados_sensibilidade["a2"] = resultados_a2
    
    # Análise do coeficiente b₁ = 1.71
    println("Analisando sensibilidade do coeficiente b₁ = 1.71...")
    resultados_b1 = SimulationResult[]
    
    for var in variacoes
        fator = 1 + var/100
        b1_novo = b1_orig * fator
        params_mod = ReactorParams(a1_orig, a2_orig, b1_novo)
        
        t, g, f = runge_kutta_4(g0, f0, dt_base, t_max, params_mod)
        label = var == 0 ? "Original (b₁=1.71)" : @sprintf("b₁=%.2f (%+d%%)", b1_novo, var)
        resultado = SimulationResult(t, g, f, params_mod, dt_base, label)
        push!(resultados_b1, resultado)
    end
    resultados_sensibilidade["b1"] = resultados_b1
    
    return resultados_sensibilidade
end

# Função para plotar análise de sensibilidade
function plotar_sensibilidade(resultados_sensibilidade)
    """
    Gera gráficos da análise de sensibilidade
    """
    println("Gerando gráficos de análise de sensibilidade...")
    
    cores = [:red, :orange, :black, :blue, :green]
    
    # Gráficos para cada coeficiente
    plots_g = []
    plots_f = []
    
    coef_nomes = ["a₁ (coef. linear)", "a₂ (coef. quadrático)", "b₁ (coef. produção)"]
    coef_keys = ["a1", "a2", "b1"]
    
    for (idx, (key, nome)) in enumerate(zip(coef_keys, coef_nomes))
        resultados = resultados_sensibilidade[key]
        
        # Gráfico g(t)
        p_g = plot(title="Sensibilidade de g(t) - $nome", 
                   xlabel="Tempo (t)", ylabel="g(t)",
                   legend=:bottomright, grid=true)
        
        # Gráfico f(t)  
        p_f = plot(title="Sensibilidade de f(t) - $nome", 
                   xlabel="Tempo (t)", ylabel="f(t)",
                   legend=:bottomright, grid=true)
        
        for (i, res) in enumerate(resultados)
            cor = cores[i]
            estilo = occursin("Original", res.label) ? :solid : :dash
            largura = occursin("Original", res.label) ? 3 : 2
            
            plot!(p_g, res.t, res.g, 
                  label=res.label, color=cor, 
                  linestyle=estilo, linewidth=largura)
            
            plot!(p_f, res.t, res.f, 
                  label=res.label, color=cor, 
                  linestyle=estilo, linewidth=largura)
        end
        
        push!(plots_g, p_g)
        push!(plots_f, p_f)
        
        # Salvar gráficos individuais
        savefig(p_g, "sensibilidade_$(key)_g.png")
        savefig(p_f, "sensibilidade_$(key)_f.png")
    end
    
    # Gráficos combinados
    p_combined_g = plot(plots_g..., layout=(3,1), size=(800, 1200))
    p_combined_f = plot(plots_f..., layout=(3,1), size=(800, 1200))
    
    savefig(p_combined_g, "sensibilidade_completa_g.png")
    savefig(p_combined_f, "sensibilidade_completa_f.png")
    
    println("Gráficos de sensibilidade salvos:")
    for key in coef_keys
        println("• sensibilidade_$(key)_g.png")
        println("• sensibilidade_$(key)_f.png")
    end
    println("• sensibilidade_completa_g.png")
    println("• sensibilidade_completa_f.png")
    
    return plots_g, plots_f
end

# Função para análise quantitativa da sensibilidade
function relatorio_sensibilidade(resultados_sensibilidade)
    """
    Gera relatório detalhado da análise de sensibilidade
    """
    open("analise_sensibilidade_completa.txt", "w") do file
        write(file, "ANÁLISE DE SENSIBILIDADE DOS COEFICIENTES\n")
        write(file, "="^70 * "\n\n")
        
        write(file, "VALORES ORIGINAIS:\n")
        write(file, "• a₁ = 13.1 (coeficiente linear)\n")
        write(file, "• a₂ = 13.94 (coeficiente quadrático)\n") 
        write(file, "• b₁ = 1.71 (coeficiente de produção)\n\n")
        
        # Valor original de referência
        res_original = findfirst(r -> occursin("Original", r.label), 
                                resultados_sensibilidade["a1"])
        g_ref = resultados_sensibilidade["a1"][res_original].g[end]
        f_ref = resultados_sensibilidade["a1"][res_original].f[end]
        
        write(file, @sprintf("VALORES DE REFERÊNCIA (t=1):\n"))
        write(file, @sprintf("• g(1) = %.8f\n", g_ref))
        write(file, @sprintf("• f(1) = %.8f\n\n", f_ref))
        
        # Análise para cada coeficiente
        coef_nomes = Dict("a1" => "a₁ (coeficiente linear)", 
                         "a2" => "a₂ (coeficiente quadrático)", 
                         "b1" => "b₁ (coeficiente de produção)")
        
        for (key, nome) in coef_nomes
            write(file, "SENSIBILIDADE DO $nome:\n")
            write(file, "-"^50 * "\n")
            write(file, @sprintf("%-15s %-12s %-12s %-10s %-10s\n", 
                  "Variação", "g(1)", "f(1)", "Δg(%)", "Δf(%)"))
            write(file, "-"^65 * "\n")
            
            for res in resultados_sensibilidade[key]
                if !occursin("Original", res.label)
                    delta_g_pct = ((res.g[end] - g_ref) / g_ref) * 100
                    delta_f_pct = ((res.f[end] - f_ref) / f_ref) * 100
                    
                    # Extrair variação do label
                    if occursin("(+", res.label)
                        variacao = split(split(res.label, "(+")[2], "%")[1] * "%"
                        variacao = "+" * variacao
                    elseif occursin("(-", res.label)  
                        variacao = split(split(res.label, "(-")[2], "%")[1] * "%"
                        variacao = "-" * variacao
                    else
                        variacao = "Original"
                    end
                    
                    write(file, @sprintf("%-15s %-12.8f %-12.8f %-10.2f %-10.2f\n", 
                          variacao, res.g[end], res.f[end], delta_g_pct, delta_f_pct))
                end
            end
            write(file, "\n")
        end
        
        write(file, "INTERPRETAÇÃO DOS RESULTADOS:\n")
        write(file, "• Δg(%) e Δf(%): variação percentual em relação ao caso original\n")
        write(file, "• Valores positivos indicam aumento da concentração\n")  
        write(file, "• Valores negativos indicam diminuição da concentração\n")
    end
    
    println("Relatório de sensibilidade salvo: analise_sensibilidade_completa.txt")
end

# FUNÇÃO PRINCIPAL
function main()
    """
    Função principal que executa toda a simulação conforme especificação
    """
    println("Iniciando simulação completa do reator de batelada...")
    println("Linguagem: Julia")
    println("Data: $(Dates.now())")
    
    # PARÂMETROS DO PROBLEMA
    g0 = 0.03     # g(0) = 0.03
    f0 = 0.0      # f(0) = 0.0  
    t_max = 1.0   # 0 ≤ t ≤ 1
    dt_inicial = 0.05  # Δt inicial = 0.05
    
    # 1. ANÁLISE DE CONVERGÊNCIA
    println("\n1️ EXECUTANDO ANÁLISE DE CONVERGÊNCIA...")
    resultados_convergencia = analise_convergencia_completa(g0, f0, t_max)
    
    # 2. PLOTAR RESULTADOS PRINCIPAIS 
    println("\n2️ GERANDO GRÁFICOS DOS RESULTADOS...")
    plotar_resultados_principais(resultados_convergencia)
    
    # 3. ANÁLISE DE SENSIBILIDADE 
    println("\n3️ EXECUTANDO ANÁLISE DE SENSIBILIDADE...")
    resultados_sensibilidade = analise_sensibilidade_coeficientes(g0, f0, t_max, 0.01)
    
    # 4. PLOTAR SENSIBILIDADE
    println("\n4️ GERANDO GRÁFICOS DE SENSIBILIDADE...")
    plotar_sensibilidade(resultados_sensibilidade)
    
    # 5. RELATÓRIO DE SENSIBILIDADE
    println("\n5️GERANDO RELATÓRIO DE SENSIBILIDADE...")
    relatorio_sensibilidade(resultados_sensibilidade)
    
    # 6. RESUMO FINAL
    println("\n" * "="^80)
    println("SIMULAÇÃO CONCLUÍDA COM SUCESSO!")
    println("="^80)
    
    println("\n ARQUIVOS GERADOS:")
    println("\n TABELAS DE RESULTADOS:")
    println("• convergencia_numerica.txt - Análise de convergência")
    println("• resultados_dt_0.0500.txt - Resultados com Δt=0.05 (principal)")
    println("• resultados_dt_*.txt - Resultados para outros valores de Δt")
    println("• analise_sensibilidade_completa.txt - Análise detalhada")
    
    println("\n GRÁFICOS PRINCIPAIS:")
    println("• evolucao_g_tempo.png - Evolução de g(t)")
    println("• evolucao_f_tempo.png - Evolução de f(t)")
    println("• plano_fase.png - Plano de fase g vs f")
    println("• resultados_principais.png - Gráficos combinados")
    
    println("\nGRÁFICOS DE SENSIBILIDADE:")
    println("• sensibilidade_a1_g.png, sensibilidade_a1_f.png")
    println("• sensibilidade_a2_g.png, sensibilidade_a2_f.png")  
    println("• sensibilidade_b1_g.png, sensibilidade_b1_f.png")
    println("• sensibilidade_completa_g.png, sensibilidade_completa_f.png")
    
    # Resultados finais com Δt = 0.05
    resultado_principal = resultados_convergencia[1]  # Δt = 0.05
    println("\n RESULTADOS FINAIS (Δt = 0.05):")
    println(@sprintf("• g(1) = %.8f", resultado_principal.g[end]))
    println(@sprintf("• f(1) = %.8f", resultado_principal.f[end]))

end

main()