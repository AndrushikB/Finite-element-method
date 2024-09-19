import numpy as np
import matplotlib.pyplot as plt


def gauss(numof_gauss_point):
    
    if numof_gauss_point == 1:
        integr = np.array([[0, 2]])
    elif numof_gauss_point == 2:
        integr = np.array([[-1./np.sqrt(3.), 1.], [1./np.sqrt(3.), 1.]])
    return integr

def find_elem(i):
#Функция возвращает координаты узлов и вектор g
    l = 0
    place = np.zeros((local_numofnodes, node_dof))
    d = np.zeros((local_numofnodes * node_dof, 1))
    for m in range(local_numofnodes):
        for n in range(node_dof):
            place[m, n] = coord[int(connect[i, m]), n]
            d[l] = dof[int(connect[i, m]), n]
            l += 1
    return place, d

def shape(integr, s, t):
#Функция возвращает вектор функций формы и их производных по xi и eta
    xi = integr[s, 0]
    eta = integr[t, 0]
    shape_fun = 0.25 * np.array([[(1. - xi) * (1. - eta)],
                               [(1. + xi) * (1. - eta)],
                               [(1. + xi) * (1. + eta)], 
                               [(1. - xi) * (1. + eta)]])
    
    deriv = 0.25 * np.array([[-(1. - eta), (1. - eta), (1. + eta), -(1. + eta)], 
                            [-(1. - xi), -(1. + xi), (1. + xi), (1. - xi)]])
    return deriv, shape_fun

def form_B(derivative, local_numofnodes, element_dof):
#Функция формирует матрицу B из производных функций формы в ГСК
    B = np.zeros((3, element_dof))
    for r in range(local_numofnodes):
        p = 2 * r + 1
        s = p - 1
        x = derivative[0, r]
        B[0, s] = x
        B[2, p] = x
        y = derivative[1, r]
        B[1, p] = y
        B[2, s] = y
    return B

def form_K(K, ke, d):
#Функция сборки глобальной матрицы жесткости
    for i in range(element_dof):
        if d[i] != 0:
            for j in range(element_dof):
                if d[j] != 0:
                    K[int(d[i, 0]) - 1, int(d[j, 0]) - 1] = K[int(d[i, 0]) - 1, int(d[j, 0]) - 1] + ke[i, j]
    return K

def average_stresses(sigma):
#Функция возвращает узловые усредненные значения напряжений    
    sigmax = np.zeros(((gridx + 1) * (gridy + 1), 1))
    sigmay = np.zeros(((gridx + 1) * (gridy + 1), 1))
    sigmaxy = np.zeros(((gridx + 1) * (gridy + 1), 1))
    for m in range(number_of_nodes):
        itr = 0 #встречаемость узла
        sx = 0
        sy = 0
        sxy = 0
        for el in range(elem):
            for node in range(local_numofnodes):
                if connect[el, node] == m:
                    itr = itr + 1
                    sx = sx + sigma[el, 0]
                    sy = sy + sigma[el, 1]
                    sxy = sxy + sigma[el, 2]
        sigmax[m, 0] = sx / itr
        sigmay[m, 0] = sy / itr
        sigmaxy[m, 0] = sxy / itr
    return sigmax, sigmay, sigmaxy

length = 60. #Длина пластины, мм
width = 20. #Ширина пластины, мм
gridx = 30 #Количество столбцов по Ox
gridy = 10 #Количество столбцов по Oy
sizex = length / gridx #Размер столбца по Ox
sizey = width / gridy #Размер столбца по Oy
OX = 0 #Начало координат
OY = width/2
local_numofnodes = 4 #Число узлов одного элемента
node_dof = 2 #Степени свободы узла
element_dof = local_numofnodes * node_dof #Степени свободы элемента

E = 200000. #Модуль Юнга
poiss = 0.3 #Коэффициент Пуассона
thick = 5. #Толщина пластины


#Построение сетки
coord = np.zeros(((gridx + 1) * (gridy + 1), 2))
connect = np.zeros((gridx * gridy, local_numofnodes))
elem = 0
number_of_nodes = 0
for i in range(gridx):
    for j in range(gridy):
        node0 = j + i*(gridy + 1)
        coord[node0, :] = [i*sizex - OX, j*sizey - OY]
        node1 = j + (i + 1)*(gridy + 1)
        coord[node1, :] = [(i + 1)*sizex - OX, j*sizey - OY]
        node2 = node1 + 1
        coord[node2, :] = [(i + 1)*sizex - OX, (j + 1)*sizey - OY]
        node3 = node0 + 1
        coord[node3, :] = [i*sizex - OX, (j + 1)*sizey - OY]
        connect[elem, :] = [node0, node1, node2, node3] #Связь узлов по КЭ
        number_of_nodes = node2 + 1 #Номер правого верхнего узла(их количество)
        elem += 1


#Матрица упругих постоянных D
D = E / (1. - poiss*poiss) * np.array([[1., poiss, 0.], 
                                      [poiss, 1., 0.], 
                                      [0., 0., (1.-poiss) / 2]])

dof = np.ones((number_of_nodes, node_dof)) # Матрица степеней свободы по узлам
#Задаем граничные условия(правая грань пластины)
for i in range(number_of_nodes):
    if coord[i, 0] == length:
        dof[i] = [0, 0]


k = 0
for i in range(number_of_nodes):
    for j in range(node_dof):
        if dof[i, j] != 0:
            k = k + 1
            dof[i, j] = k
            


#Задаем нагрузки
P = 1000 #Узловая сила
force = np.zeros((number_of_nodes, 2))
for i in range(number_of_nodes):
    if coord[i, 0] == 0 and coord[i, 1] == 0: #Сила в узле (0,0)
        force[i, :] = [0, -P]


forcevect = np.zeros((k, 1))
for i in range(number_of_nodes):
    if dof[i, 0] != 0:
        forcevect[int(dof[i, 0]) - 1] = force[i, 0]
    if dof[i, 1] != 0:
        forcevect[int(dof[i, 1]) - 1] = force[i, 1]
        
     
numof_gauss_point = 2
integr = gauss(numof_gauss_point) #Гауссовы точки и весовые коэффициенты


#Вычисление глобальной матрицы жесткости
K = np.zeros((k, k))
for i in range(elem):
    place, d = find_elem(i)
    ke = np.zeros((element_dof, element_dof)) #Матрица жесткости элемента 8*8
    #Весовые коэффициенты, квадратуры Г.-Л.
    for s in range(numof_gauss_point):
        wi = integr[s, 1]
        for t in range(numof_gauss_point):
            wj = integr[t, 1]
            deriv, shape_fun = shape(integr, s, t) #Произв. функций формы в ЛСК
            jacobi = np.dot(deriv, place)
            detj = np.linalg.det(jacobi)
            inv_jacobi = np.linalg.inv(jacobi) #Обратная матрица Якоби
            derivative = np.dot(inv_jacobi, deriv) #Произв. функций формы в ГСК
            B = form_B(derivative, local_numofnodes, element_dof)
            ke = ke + detj * thick * wi * wj * np.dot(np.dot(B.T, D), B) #Вычисление матрицы жесткости КЭ
            #print(ke)
            #time.sleep(2)
    K = form_K(K, ke, d)


delta = np.linalg.solve(K, forcevect)


displacements = np.ones(((gridx + 1) * (gridy + 1), 2))

for i in range(number_of_nodes):
    if dof[i, 0] == 0:
        xmove = 0
    else:
        xmove = delta[int(dof[i, 0]) - 1, 0]

    if dof[i, 1] == 0:
        ymove = 0
    else:
        ymove = delta[int(dof[i, 1]) - 1, 0]

    displacements[i, :] = [xmove, ymove]

numof_gauss_point = 1
integr = gauss(numof_gauss_point)
sigma = np.zeros((gridx * gridy, 3))
for i in range(elem):
    place, d = find_elem(i)
    vectdisp = np.zeros((element_dof, 1))
    for m in range(element_dof):
        if d[m] == 0:
            vectdisp[m] = 0
        else:
            vectdisp[m] = delta[int(d[m, 0]) - 1]

    for s in range(numof_gauss_point):
        wi = integr[s, 1]
        for t in range(numof_gauss_point):
            wj = integr[t, 1]
            deriv, shape_fun = shape(integr, s, t)
            jacobi = np.dot(deriv, place)
            detj = np.linalg.det(jacobi)
            inv_jacobi = np.linalg.inv(jacobi)
            derivative = np.dot(inv_jacobi, deriv)
            B = form_B(derivative, local_numofnodes, element_dof)
            epsilon1 = np.dot(B, vectdisp) #деформация
            sigma1 = np.dot(D, epsilon1) #напряжения

    sigma[i, :] = sigma1.T

sigmax, sigmay, sigmaxy = average_stresses(sigma)
Eg = np.dot(np.dot(delta.T, K), delta / 2)
v = displacements[:, 1]
u = displacements[:, 0]
cmin = min(sigmay)
cmax = max(sigmay)


#Смещение по OX
X, Y = np.meshgrid(np.linspace(0, length, gridx+1), np.linspace(-width/2, width/2, gridy+1))
Z = u.reshape(gridy+1, gridx+1, order = 'F')
levels = np.linspace(np.min(Z), np.max(Z), 20)
c1 = plt.contour(X, Y, Z, levels = levels)


c1.clabel()
plt.grid(which = 'major')
plt.title("Displacement Ox")
plt.xlabel("Ox", fontstyle = 'italic')
plt.ylabel("Oy", fontstyle = 'italic')
#matplotlib.colors.Liste

plt.show(c1)

#Смещение по OY
X, Y = np.meshgrid(np.linspace(0, length, gridx+1), np.linspace(-width/2, width/2, gridy+1))
Z = v.reshape(gridy+1, gridx+1, order = 'F')
levels = np.linspace(np.min(Z), np.max(Z), 15)
c2 = plt.contour(X, Y, Z, levels = levels)

c2.clabel()
plt.grid(which = 'major')
plt.title("Displacement Oy")
plt.xlabel("Ox", fontstyle = 'italic')
plt.ylabel("Oy", fontstyle = 'italic')

plt.show(c2)


      