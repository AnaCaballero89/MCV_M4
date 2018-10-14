"# MCV_M2" 
TASKS TO DO:
==============================================
======================= MANDATORY ==============
==============================================
Mandatory means if there any point that it is not done, then the weekly task will FAIL. Implement the importing gradients method from the Patrick’s paper. 

* Read the Patrick Perez’s paper. 
* Complete start.m file 
-> he sustituido sol_ por G5 en todo, y he rellenado el fichero de start.m
-> hay un if para elegir la foto

* Modify the G??? Laplace Equation Axb.m file and create a G??? Poisson Equation Axb.m adding ONLY 4 lines. 3 of this 4 lines should be 
 - if (isfield(param, ’driving’)) 
 - else 
 - end
-> hecho pero no son 4, son 5. En el pixek que cambia:
            if (isfield(param, 'driving'))
                b(p) = param.driving(i-1,j-1);
            else
                b(p) = 0;
            end 
* Test with your own images
-> He añadido 2 imágenes y dos fondos(toinsertX y backgroundX) pero esta sin probar.
-> la primera es un pato en una cascada y la otra una ardilla en un jardin.
-> también podemos jugar con nuestras caras y las dos que nos dan (estoy en ello, no me funciona bien)

* The objective is that if param.driving exists, do the Poisson editing, and if it doesn’t, do the Inpainting by Laplace’s equation. 
* Deliverable of mandatory. Thu. Oct. 18. 18h 
* WARNING: Be careful with the sign on the discretization of the Laplacian operator!!!!!!

============================================== 
======================= OPTIONAL ==============
============================================== 
1. Test with your own real images. Each successfully edited image will sum up 0.1 points. 
* Cuando funcione la ultima del opcional tendremos 11 combinaciones, incluida la que nos dan.

2. +1 Points: Implement Seamless cloning with mixing gradients. 
3. Try different numerical schemes 
3.1 +5 points: Gradient Descent and compare it with cyclic schemes for evolution PDEs 
3.2 +5 points: Gauss-Seidel and ω-relaxation and compare them with cyclic schemes for stationary PDEs 

3.3 +10 points: Multigrid. 

4. Feel free to experiment
