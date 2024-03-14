---
title: "Análise Fatorial - Modelo Via Matriz de Correlação no R"
author: "Gabriel Carvalho, Gustavo Henrique, Mayara Farias, Tamires Domingos."
date: "12/03/2024"
bibliography: references.bib
output:
  html_document:
    df_print: paged
  word_document: default
  pdf_document:
    keep_tex: true
header-includes:
- \usepackage{fancyhdr}
- \pagestyle{fancy}
- \fancyhead{}
- \fancyfoot{}
- \fancyfoot[R]{\thepage}
---

```{r}
# Carregar o conjunto de dados
db <- read.csv("https://raw.githubusercontent.com/gabrieldacarvalho/analise_fatorial/main/RTE_CEREAL.csv")

# Visualizar as primeiras linhas do conjunto de dados
print(head(db))

# Selecionar as variáveis para análise (excluindo as variáveis de identificação)
x <- as.matrix(db[, 3:27])

# Obter o número de observações
n <- nrow(x)

# Calcular a média de cada variável
matriz_media <- matrix(1, nrow = 1, ncol = n)
matriz_media <- (matriz_media %*% x) / n 

# Centralizar os dados subtraindo a média de cada variável
x_centralizado <- sweep(x, 2, matriz_media, "-")

# Padronizar os dados dividindo pelo desvio padrão de cada variável
z <- sweep(x_centralizado, 2, apply(x, 2, sd), "/")

# Calcular a matriz de correlação amostral
matriz_r <- (t(z) %*% z) / (n - 1)

```

## Estimação do número de fatores m:

Para estimar o número ideal de fatores, procedemos da seguinte maneira:

(1) Estimar a matriz de correlação teorica \(P_{p \times p}\) a partir da matriz de correlação amostral \(R_{p\times p}\).

(2) Obter os autovalores da matriz \(R_{p\times p}\)

(3) Ordenar os autovalores em ordem crescente e pegar a quantidade de autovalores
maior ou igual a 1 para ser a quantidade de fatores m

```{r}
# Visualizar os autovalores e autovetores
eigen_vals <- eigen(matriz_r)$values
eigen_vecs <- eigen(matriz_r)$vectors

print("Autovalores:")
print(eigen_vals)
print("Autovetores:")
print(eigen_vecs)

# Selecionar os m maiores autovalores e autovetores correspondentes
m <- sum(eigen_vals > 1)
autovalores <- eigen_vals[1:m]
autovetores <- eigen_vecs[, 1:m]
  
``` 

### Estimação para matriz L e ψ

```{r}
matriz_l <- sqrt(autovalores) * autovetores
matriz_psi <- matriz_r - matriz_l %*% t(matriz_l)
matriz_r_hat <- matriz_l %*% t(matriz_l) + matriz_psi

# Visualizar as matrizes estimadas
print("Matriz L:")
print(matriz_l)
print("Matriz ψ:")
print(matriz_psi)
print("Matriz estimada R:")
print(matriz_r_hat)

# Calcular a matriz residual
matriz_residual <- matriz_r - matriz_r_hat

# Visualizar a matriz residual e a soma dos resíduos
print("Matriz Residual:")
print(matriz_residual)
print("Soma dos Resíduos:")
print(sum(matriz_residual))
```

A estimativa derivada do modelo ortogonal demonstrou ser altamente precisa, resultando em uma matriz residual significativamente reduzida.   


### Estimação dos Escores dos Fatores F

Para calcular o escores dos fatores vamos utilizar o método dos mínimos quadrados
ponderados
\[
F = (L'\ ψ^{-1}\ L)^{-1}\ L'\ ψ^{-1}\ Z
\]

```{r}
# Calcular a matriz F
matriz_f <- solve(t(matriz_l) %*% solve(matriz_psi) %*% matriz_l) %*% t(matriz_l) %*% solve(matriz_psi) %*% t(z)

# Imprimir a soma das linhas dos escores dos fatores
print(matrix(rowSums(matriz_f)))

# Imprimir a matriz resultante da multiplicação entre a matriz L e a matriz F
print(t(matriz_l %*% matriz_f)[1:5,])
print(z[1:5,])
```


# Referências
Silva, Wilton Bernardino da. 2021. “Análise Fatorial (Modelo Ortogonal).” YouTube. 2021. https://www.youtube.com/watch?v=Jk-QGfBkifA.
