import numpy
import sequential_lib as seq

size = int(input("Matrix size: "))

A = numpy.random.rand(size,size)
B = numpy.random.rand(size,size)

expected = numpy.matmul(A, B)

result, status = seq.mult_seq_explicit(A, B)
print("Result:")
print(result)
print("Expected:")
print(expected)

A = numpy.array(numpy.random.rand(size,size), order='F')
X = numpy.array(numpy.random.rand(size), order='F')

print("Input:")
print(A)
print(X)

result = seq.gauss_seq(A, X, size - 1)

print("Result: ")
print(numpy.array(A))
print(" ")
print(numpy.array(X))
