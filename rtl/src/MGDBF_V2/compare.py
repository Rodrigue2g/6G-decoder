## read xhat_out.txt and y_out.txt and compare them

import numpy as np

y= np.loadtxt('y_out.txt') 
xhat = np.logical_not(np.loadtxt('xhat_out.txt'))
xhat_matlab = np.loadtxt('xhat_out_MATLAB.txt')
x_matlab = np.loadtxt('x_out_MATLAB.txt')


errors_matlab = 0
errors_real = 0
errors_C = 0
errors_HD=0
for i in range(len(xhat)-640):

    if xhat[i] != (y[i]<0):
        errors_HD += 1

    if x_matlab[i] != xhat_matlab[i]:
        errors_C += 1

    if xhat[i] != x_matlab[i]:
        errors_real += 1

    if xhat_matlab[i] != xhat[i]:
        errors_matlab += 1
        print(int(xhat[i]), int(xhat_matlab[i]), int(x_matlab[i]), int(y[i]<0), y[i], "                  ERROR!")
    else:
       print(int(xhat[i]), int(xhat_matlab[i]), int(x_matlab[i]), int(y[i]<0), y[i])

print(np.max(np.loadtxt('y_out.txt')))
print("Errors vivado - decoded matlab: ", errors_matlab)
print("Errors vivado - sent matlab: ", errors_real)
print("Errors decoded matlab - sent matlab: ", errors_C)
print("Errors Hard Decision: ", errors_HD)
