# A interface (?)/ data structure 

* shared library libdengo.so
* data structure that can be passed around by the simulations
* pointers to where the field data (HI, H2I ...) is
* modules to calculate temperature, pressure based on the field data
* **solver** modules that takes a field_data object, dt, and return 

\begin{equation}
y = f(x)
\end{equation}