using ITensors

# Example: create a simple ITensor and print it
i = Index(2, "i")
T = randomITensor(i)
println("Random ITensor:")
println(T)