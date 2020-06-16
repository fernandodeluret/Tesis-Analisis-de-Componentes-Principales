![](media/image1.jpg){width="4.365in" height="0.6866666666666666in"}

> Carrera: Licenciatura en Economía

ÍNDICE DE ACTIVIDAD ECONÓMICA DE MENDOZA, POR ANÁLISIS DE COMPONENTES PRINCIPALES 
==================================================================================

Trabajo de Investigación

POR

Fernando Deluret

Profesor Tutor

Mónica Calderón

M e n d o z a - a ñ o 2 0 1 9

INDICE 
=======

  --------------------------------------------------------------------------------------
  Introducción .………………………………………………………………………………………………….……………………........              3
  --------------------------------------------------------------------------------- ----
  CAPÍTULO I – CONCEPTOS PREVIOS …………………………………………………………………….……………….........         4

  1.  Vectores y espacios vectoriales ………………………………………………………………….………….               4
                                                                                    

  1.  Ortogonalidad ……………………..……………………………………………………………….……………...                     5
                                                                                    

  1.  Proyecciones ortogonales ……………………………………………………………………………………..                   5
                                                                                    

  1.  Complemento orthogonal …………………………………………………………….……………………….                     8
                                                                                    

  CAPÍTULO II – ANÁLISIS DE COMPONENTES PRINCIPALES …………………………………………………….........   10

  CAPÍTULO III – ANÁLISIS FACTORIAL …………………………………………………………………………………….........       19

  CAPÍTULO IV – ROTACIÓN DE COMPONENTES …………………………………………………………………….........         23

  CAPÍTULO V – CASO PRÁCTICO …………………………………………………………………………………………….........           25

  Conclusiones ………………………………………………………………………………………………….……………………........               32

  Bibliografía consultada ….……………………………………………………………………………….……………………........         34
  --------------------------------------------------------------------------------------

INTRODUCCIÓN

El análisis de componentes principales (ACP) es un método de reducción
de dimensionalidad lineal de los datos. El objetivo es poder reducir el
número de variables necesario para representar un conjunto de datos.
Para esto lo que se hace es realizar una proyección ortogonal de los
datos originales en un subespacio vectorial de menor dimensión, y que
por ende pueda ser expresado como un conjunto menor de variables.

El presente trabajo se estructurará de la siguiente manera, en el
capítulo 1 se analizarán una serie de conceptos previos necesarios tanto
para el desarrollo matemático del ACP como para su interpretación
práctica. En el capítulo 2 se realizará el desarrollo a partir del cual
surge la fórmula resultante del ACP. En el capítulo 3 se analizará el
tema desde la perspectiva del análisis factorial que nos permite añadir
el componente de interpretabilidad de los componentes. En el capítulo 4
se analizará la rotación de componentes, que es la herramienta que nos
permitirá relacionar las variables originales con los componentes o
factores interpretando estas últimas como variables subyacentes
observadas a través de los datos iniciales. En el capítulo 5 se
trabajará con un caso numérico, el cual tendrá por objetivo tratar de
construir un índice de actividad mensual para la provincia de Mendoza.
Y, por último, se analizarán las consecuencias prácticas del ACP y se
sacarán algunas conclusiones.

**CAPÍTULO I**

**CONCEPTOS PREVIOS**

1.  Vectores y espacios vectoriales

Trabajando en $\mathbb{R}^{2}$ y dados 2 vectores
$\overline{\text{AB}}\ $ y $\overline{\text{AC}}$ como se muestra en el
Gráfico 1, podemos definir el módulo de un vector, la distancia entre
vectores y el ángulo entre vectores de la siguiente manera:

*Gráfico 1*

> ![](media/image2.png){width="3.1791666666666667in"
> height="3.1458333333333335in"}$\text{Siendo}\text{\ \ }w_{1} = \overline{\text{AB}}\text{\ \ y\ \ }w_{2} = \overline{\text{AC}},\ \text{con}\ w_{i} = \begin{bmatrix}
> x_{i} \\
> y_{i} \\
> \end{bmatrix}$
>
> Podemos definir el producto punto entre 2 vectores como:
>
> ${w_{1}}^{T}\text{.\ }w_{2} = \ x_{1}\text{.x}_{2} +$ $y_{1}.y_{2}$
>
> El *módulo* de $w_{1}$ será $||w_{1}|| = \sqrt{w^{T}\text{.w}}$ y esto
> es igual al segmento $\overline{\text{AB}}$.
>
> La *distancia* entre $w_{1}\text{\ y\ }w_{2}$ será el modulo del
> vector diferencia entre ambos $||w_{1} - w_{2}||$ lo cual va a ser
> igual al segmento $\overline{\text{BC}}$.
>
> El *ángulo* $\alpha$ entre estos 2 vectores va a ser tal que:
> $\cos{(\alpha)} = \frac{{w_{1}}^{T}\text{.\ }w_{2}}{||w_{1}||.||w_{2}||}$

Un espacio vectorial es un conjunto de vectores que cumplen las
siguientes propiedades:

-   Son cerrados respecto de una operación (por ejemplo la suma): si
    $w_{1},w_{2} \in V\  \Rightarrow \ w_{1} + w_{2} \in V$

-   Asociatividad:
    $w_{1} + {(w}_{2} + w_{3}) = {(w}_{1} + w_{2}) + w_{3}$

-   Elemento neutral (*e*):
    $\exists e\  \in V\ \ \ \text{tal\ que}\text{\ \ \ \ \ \ \ \ }w_{1} + e = {e + w}_{1} = w_{1}$

-   Elemento inverso:
    $\exists{w_{1}}^{- 1}\  \in V\ \ \text{tal\ que}\text{\ \ \ \ \ \ \ \ }w_{1} + {w_{1}}^{- 1} = e$

-   Multiplicación por un escalar:
    $w_{1} \in V\  \Rightarrow \ \lambda.w_{1} \in V$

Un subespacio vectorial es entonces un subconjunto de ese espacio tal
que si aplicamos cualquiera de las operaciones arriba mencionadas entre
vectores del subespacio el vector resultante seguirá perteneciendo al
mismo. Una consecuencia de esto es que el subespacio debe tener si o si
el elemento neutro, por ejemplo, las rectas que pueden ser subespacio de
$\mathbb{R}^{2}$ deben pasar si o si por el origen de coordenadas.

Dado un espacio vectorial V, un *conjunto generador* para este espacio
es un conjunto de vectores que combinados linealmente permiten construir
cualquier vector de dicho espacio. Si adicionalmente el conjunto
generador cumple las siguientes propiedades se dice que este conjunto es
una *base* de ese espacio vectorial y a cada uno de los vectores de la
base se le llama *vector base:*

-   Mínimo:
    $A = \ \left\lbrack w_{1},w_{2},\ \ldots\ ,w_{k} \right\rbrack\ \ \ \ A \subseteq V$,
    será mínimo si no existe otro conjunto generador con menos de k
    vectores.

-   Los vectores de A deben ser linealmente independientes

Una propiedad de una base es que 2 combinaciones lineales diferentes de
sus vectores base darán como resultado 2 vectores diferentes entre sí
siempre.

Una *base estándar* de un espacio $\mathbb{R}^{n}$, es una base tal que
la matriz formada por sus vectores base es la matriz identidad de
dimensión $n \times n$. Por ejemplo, la base estándar para
$\mathbb{R}^{3}$ es:

$$A = \ \left\{ \underset{\begin{matrix}
\text{vector} \\
\text{base} \\
\end{matrix}}{}\ ,\ \begin{bmatrix}
0 \\
1 \\
0 \\
\end{bmatrix},\ \begin{bmatrix}
0 \\
0 \\
1 \\
\end{bmatrix} \right\}$$

Puede haber otras bases que no sean la estándar para un espacio
vectorial, pero todas van a tener la misma cantidad de vectores.

1.  Ortogonalidad

Se dice que dos vectores son ortogonales entre sí, si el ángulo que
forman es de 90°. Por lo que, dada la fórmula del ángulo entre vectores:

> $\cos{(90^{o})} = \frac{{w_{1}}^{T}\text{.\ }w_{2}}{||w_{1}||.||w_{2}||} = 0$

Para que esto se cumpla el producto punto entre ambos vectores debe ser
igual a cero. Una forma de interpretar la ortogonalidad es que dos
vectores ortogonales son “lo más diferentes posible” entre sí.

Si los vectores de una base son ortogonales entre si
${w_{i}}^{T}\text{.\ }w_{k} = 0\ \ \ i \neq k\ \ $y todos tienen módulo
igual a 1 $\left| \left| w_{i} \right| \right| = 1$ entonces la base es
*ortonormal*.

1.  Proyecciones ortogonales

Dado el subespacio vectorial $U \subset \mathbb{R}^{2}$ y el vector
$w \in \mathbb{R}^{2}$, el cual puede ser representado como una
combinación lineal de los vectores base de $\mathbb{R}^{2}$. Supongamos
que queremos representar $w$ en el subespacio U, vamos a querer que el
nuevo vector $w' \in U$ sea lo más parecido al vector original. Lo cual
es equivalente a minimizar el segmento que representa el vector
diferencia $(w - w')$, el cual en el Gráfico 2 está representado por las
distintas líneas punteadas para cada caso.

Siendo $b_{1} = \ \begin{bmatrix}
b_{11} \\
b_{12} \\
\end{bmatrix} \in U$ una base del subespacio U, el nuevo vector $w'$ va
a poder ser expresado como una combinación lineal de dicha base
$w^{'} = \ \beta.b_{1}$, donde $\beta$ es un escalar y representa las
coordenadas de $w'$ en U.

La pendiente de U nos la va a dar la relación entre las componentes
$x_{1}\ e\ y_{1}$ del vector base $b_{1}$.

![](media/image3.png){width="3.6023972003499565in"
height="3.384615048118985in"}*Gráfico 2*

![](media/image4.png){width="2.828472222222222in"
height="2.504166666666667in"}

Como se puede observar en el grafico anterior de acuerdo al $\beta$ que
elijamos va a ser el largo del vector $w'$ resultante. Y por ende, el
largo del vector diferencia $(w - w')$ va a depender también de $\beta$.
Existirá un $\beta^{*}$ que minimice esa distancia, el cual va a ser tal
que haga que el vector diferencia y el vector base $b_{1}$ sean
ortogonales ($\alpha = 90^{o}$). Esto como se vio en el apartado
anterior es equivalente a decir que el producto punto entre ambos
vectores es igual a cero, matemáticamente se puede expresar:

> ${b_{1}}^{T}\text{.\ }\left( w - w^{'} \right) = 0$
>
> $\Rightarrow \left( {b_{1}}^{T}.w) - ({b_{1}}^{T}.w^{'} \right) = 0$
> Por propiedad de bilineariedad del producto punto
>
> $\Rightarrow {{(b}_{1}}^{T}.w) - ({b_{1}}^{T}\text{.β.}b_{1}) = 0$
> Usando la definición de w’

$$\Rightarrow {{(b}_{1}}^{T}.w) - \beta(\underset{{||b_{1}||}^{2}}{}) = 0$$

$$\Rightarrow \beta = \frac{{b_{1}}^{T}\text{.w}}{{||b_{1}||}^{2}}$$

$$\Rightarrow w^{'} = \beta.b_{1} = \underset{\begin{matrix}
\text{Matriz\ } \\
\text{de\ } \\
proyeccion \\
\end{matrix}}{}\text{w\ }$$

La matriz de proyección es aquella matriz (en este caso $2 \times 2$)
que multiplicada por el vector original w nos da la proyección de este
en el subespacio U (w’). Si $b_{1}$ es ortonormal
(${||b_{1}||}^{2} = 1$) la matriz de proyección va a ser simétrica:

$$b_{1}.{b_{1}}^{T} = \left( b_{1}.{b_{1}}^{T} \right)^{T}$$

Se puede ver entonces, que antes para representar al vector w
necesitábamos dos coordenadas ($x_{1}\ e\ y_{1}$) y ahora necesitamos
una sola ($\beta$).

Algo importante a destacar en este tema es que el nuevo vector w’ va a
seguir perteneciendo a $\mathbb{R}^{2}$, a pesar de que ahora solo
dependa de un parámetro ($\beta$). Esto es debido a que el vector base
de U ($b_{1}$) tiene 2 componentes y por ende w’ también. Por este
motivo no sería lo mismo por más que hagamos una proyección sobre una
recta partir desde $\mathbb{R}^{2}$, que por ejemplo desde
$\mathbb{R}^{3}$, en cuyo caso $b_{1}$ tendría 3 componentes. Una forma
de ver esto es que $b_{1}$ nos está dando información de acuerdo a los
grados de libertad de rotación del espacio vectorial inicial, por
ejemplo en $\mathbb{R}^{2}$ la recta U solo puede rotar a lo largo del
plano XY en tanto que en $\mathbb{R}^{3}$ la recta podría rotar respecto
al plano respecto a XY, YZ y XZ.

Supongamos ahora que $w \in \mathbb{R}^{D}$ es decir
$w = \left\lbrack w_{1},\ w_{2},\ldots,\ w_{D} \right\rbrack$ y
$U \subset \mathbb{R}^{D}$, pero ahora en vez de ser una recta
(dimensión igual a 1) es un subespacio M-dimensional, con $M < D$. Es
decir, tiene M vectores base
$b = \left\lbrack b_{1},\ b_{2},\ldots,\ b_{M} \right\rbrack$ donde a su
vez cada $b_{i}$ tiene dimensión $D \times 1$ ya que U es un subespacio
de $\mathbb{R}^{D}$.

Para aclarar esto último, la cantidad de vectores base $b_{i}$ que tenga
U dice si el o los vectores que se proyecten lo harán sobre una recta,
un plano etc. Y “D” lo que dice es que esa recta o plano debe dibujarse
en el espacio D-dimensional. Por ejemplo $b = \begin{bmatrix}
0 \\
0 \\
1 \\
\end{bmatrix}$ genera una recta en $\mathbb{R}^{3}$ en tanto que el
conjunto $b = \left\{ \begin{bmatrix}
0 \\
0 \\
1 \\
\end{bmatrix},\begin{bmatrix}
0 \\
1 \\
0 \\
\end{bmatrix} \right\}$ genera un plano en $\mathbb{R}^{3}$, luego sobre
ese plano o recta es sobre el que se proyectara el vector original.

En este caso multidimensional
$w^{'} = \ \sum_{i = 1}^{M}{\beta_{i}.b_{i}}$ notado matricialmente
$w^{'} = B.\beta$, donde B es la matriz que representa el conjunto de
los vectores base de U (con dimensión $D \times M$) y $\beta$ es el
vector $M \times 1$ con las M coordenadas.

La condición que se ponía antes para que la proyección fuera ortogonal
es que el vector diferencia fuera ortogonal con el vector base del nuevo
subespacio. Por ende lo que se va a exigir ahora es que $(w - w')$ sea
ortogonal con cada uno de los M vectores base. Para verlo con un
ejemplo, supóngase el caso en que se quiere proyectar un vector de
$\mathbb{R}^{3}$ en el plano XY (es decir componente z=0) el cual lo
vamos a suponer horizontal como se ilustra en el gráfico 3, en tanto que
la componente Z daría la altura. La forma de que el vector diferencia
sea lo más pequeño posible seria bajando verticalmente, es decir en
forma perpendicular al plano XY, con lo cual el vector diferencia va a
ser ortogonal a todos los vectores de XY incluidos ambos vectores base.

*Gráfico 3*

> ![D:\\Asset 8.png](media/image5.png){width="1.8374376640419947in"
> height="3.1076924759405076in"}
>
> ${b_{i}}^{T}\text{.\ }\left( w - w^{'} \right) = 0\ \ \ \ \ con\ \ i = 1,\ 2,\ \ldots,\ M$

Es decir, ahora hay un sistema de M ecuaciones simultáneas. Notado
matricialmente:

$${\underset{\text{MxD}}{}\left( \underset{Dx1}{} - \underset{Dx1}{} \right) = 0\backslash n}{\Rightarrow B^{T}.w - B^{T}.B.\beta = 0\backslash n}{\Rightarrow \beta = {\left( B^{T}\text{.B} \right)^{- 1}B}^{T}\text{.w}\backslash n}{\Rightarrow w^{'} = B.\beta = \underset{\begin{matrix}
\text{Matriz\ de\ proyeccion} \\
(DxD) \\
\end{matrix}}{}\text{.w}}$$

Donde ahora w’ y $\beta$ son las fórmulas para calcular la proyección
pero para el caso general de D dimensiones y llevar al vector a un
subespacio con solo M coordenadas ($M < D$). Una vez más w’ seguirá
siendo un vector $D \times 1$.

1.  Complemento ortogonal

Dado un espacio vectorial $V^{n}$ y un subespacio
$W \subseteq V,\ \text{con}\ W^{k}$ (siendo n y k las dimensiones de los
espacios). Entonces el complemento ortogonal de W es un subespacio
$W^{\bot}$ de dimensión (n-k), tal que $W^{\bot}$ contiene todos los
vectores de V que sean ortogonales a todos los vectores de W.

Cada vector $x \in V$ se puede *descomponer ortogonalmente* (de forma
única) de la siguiente forma:

$$x = \sum_{i = 1}^{k}{\delta_{i}b_{i}} + \sum_{i = 1}^{n - k}{\varphi_{i}{b_{i}}^{\bot}}$$

Donde $b_{i}$ y ${b_{i}}^{\bot}$ son los vectores base de los
subespacios W y $W^{\bot}$ respectivamente. Lo único que se está
diciendo es que cualquier vector $x$ se puede expresar como una suma de
un vector de W más un vector de $W^{\bot}$. Por ejemplo, supongamos que
W y $W^{\bot}$ son 2 rectas pertenecientes a $\mathbb{R}^{2}$
perpendiculares entre sí, cualquier vector de $\mathbb{R}^{2}$ puede ser
expresado como una suma de vectores como en el Gráfico 4.

*Gráfico 4*

![](media/image6.png){width="2.4653532370953632in" height="2.375in"}

**CAPÍTULO II**

**ANÁLISIS DE COMPONENTES PRINCIPALES**

Dado un conjunto de datos de D variables y N observaciones de cada
variable, podemos representar los datos como un conjunto de N vectores
pertenecientes a $\mathbb{R}^{D}$:

$${X = \left\{ x_{1},\ x_{2},\ldots,\ x_{N,} \right\}\text{\ \ \ }x_{j}\text{\ ϵ}\mathbb{\text{\ R}}}^{D}$$

Donde cada $x_{i}$ está formado por una observación de cada una de las D
variables iniciales. El análisis de componentes principales busca
encontrar una representación con una dimensión menor pero que sea lo más
parecida posible a X. Para lograr esto lo que se hace es realizar una
proyección ortogonal de los datos, ya que esto minimiza el vector
diferencia con respecto a los datos originales.

Para ver esto con un ejemplo, supongamos que como se muestra en el
gráfico 5 tenemos dos variables iniciales (D=2)
$z_{1}\text{\ y\ }z_{2}$. Podemos representar gráficamente el par de
datos como un vector en $\mathbb{R}^{2}$, entonces tendremos tantos
vectores como cantidad de observaciones haya en las series de datos.
Nuestro objetivo será entonces encontrar una representación de menor
dimensionalidad de este conjunto de datos, en este caso la única
dimensionalidad menor a dos es uno, o sea que el objetivo sería ajustar
esa serie de vectores a una recta. Existirán entonces incógnitas a
resolver, primero cómo representar cada uno de los vectores iniciales de
tal forma que haya el menor error posible, y segundo cómo elegir la
pendiente de esa recta de tal forma que se minimice el conjunto de los
errores. Es decir, minimizar la suma de los segmentos punteados del
gráfico 5. Ahora bien, recordando lo visto en proyecciones ortogonales
eso da respuesta a la primera de las incógnitas, la forma de representar
a cada vector dentro del subespacio U es proyectarlo ortogonalmente, ya
que de esta manera se asegura que se minimiza el error de proyección.
Queda entonces resolver el problema de encontrar la pendiente de la
recta U.

Una acotación metodológica a hacer en este punto, es que, si bien se
habla de la recta U, esta es en realidad un subespacio vectorial. Por lo
que no puede ser cualquier recta, debe pasar por el origen de
coordenadas. La consecuencia práctica de esto es que a la hora de
trabajar con ACP los datos deben tener media igual a cero.

*Gráfico 5*

![D:\\Asset 2.png](media/image7.png){width="5.725085301837271in"
height="5.444444444444445in"}

Para empezar recordemos que cada uno de los vectores $x_{j}$ puede ser
representado como una combinación lineal de los vectores base de
$\mathbb{R}^{D}$:

$$x_{j} = \ \sum_{i = 1}^{D}{B_{\text{ij}}\text{\ .}b_{i}}\text{\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ }\left( 1 \right)$$

Donde $b_{i}$ es cada uno de los vectores base y $B_{\text{ij}}$ es un
escalar que representa la coordenada correspondiente al vector $b_{i}$.
Cada $b_{i}$ va a tener dimensión D$x$1, ya que es un vector base de
$\mathbb{R}^{D}$. Supondremos para trabajar con ACP que las bases van a
ser ortonormales, es decir que $\sum_{i = 1}^{D}{{b_{i}}^{2} = 1}$.

$B_{\text{ij}} = {x_{j}}^{T}\text{\ .}b_{i}\text{\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ }\left( 2 \right)$

La ecuación (2) surge de la deducción de proyección ortogonal para las
coordenadas, donde $\left| {|b}_{i}| \right| = 1$ por ser las bases que
usaremos para ACP ortonormales. $B_{\text{ij}}\text{\ .}b_{i}$ puede ser
visto como la proyección ortogonal de $x_{j}$ en el subespacio
unidimensional generado por $b_{i}$ (es la proyección sobre esa recta o
eje). Una aclaración a hacer aquí es que la dimensión inicial de los
datos es D (se necesitan D coordenadas para representarlos), y lo que se
busca es reducir esa dimensión a M minimizando el error en el proceso.
La ecuación (2) muestra cual va a ser la magnitud de una de esas M
coordenadas, más específicamente la coordenada sobre el subespacio
unidimensional generado por el vector $b_{i}$.

Para el ACP se han hecho 2 supuestos importantes, uno es que los datos
deben tener media cero y el otro es que las bases van a ser
ortonormales, profundicemos un poco en este último supuesto.

El hecho que las bases sean ortonormales significa dos cosas, una ya la
dijimos y es que la norma de cada vector base es 1
($\sum_{i = 1}^{D}{{b_{i}}^{2} = 1}$), la otra es que los vectores base
van a ser ortogonales entre sí, es decir, su producto punto va a ser
igual a cero (${b_{i}}^{T}.b_{j} = 0\ \text{para}\ todo\ i \neq j\ $).
Gráficamente esto último significa que los vectores base, que pueden ser
vistos como los nuevos ejes de referencia para los datos transformados
(ya que es respecto de quienes están dadas las coordenadas) serán
perpendiculares entre sí.

Con esto en mente retomemos la fórmula obtenida para proyecciones
ortogonales M-dimensionales:

$${\widetilde{x}}_{j} = B.\beta = {B\left( B^{T}\text{.B} \right)^{- 1}B}^{T}.x_{j}$$

Donde ${\widetilde{x}}_{j}$ es la proyección ortogonal del vector
$x_{j}$ en el subespacio M-dimensional y B es la matriz formada por el
conjunto de los M vectores base de dicho subespacio
$B = \ \left\lbrack b_{1},\ b_{2},\ \ldots,\ b_{M} \right\rbrack$ donde
cada $b_{i}$ tiene dimensión $D \times 1$. Analizando la matriz
resultante de $B^{T}\text{.B}$ tendremos que:

$$B^{T}.B = \ \begin{bmatrix}
{b_{1}}^{T}.b_{1} & \ldots & {b_{1}}^{T}.b_{M} \\
 \vdots & \ddots & \vdots \\
{b_{M}}^{T}.b_{1} & \ldots & {b_{M}}^{T}.b_{M} \\
\end{bmatrix} = \ \begin{bmatrix}
1 & \ldots & 0 \\
 \vdots & \ddots & \vdots \\
0 & \ldots & 1 \\
\end{bmatrix}$$

Dado que por ser bases ortonormales
${b_{i}}^{T}.b_{j} = 0\ \text{para}\ todo\ i \neq j\ y\ {b_{i}}^{T}.b_{i} = 1$,
entonces tenemos como resultado la matriz identidad $M \times M$. Por lo
que podemos reescribir ${\widetilde{x}}_{j}$ de la siguiente manera:

$${\widetilde{x}}_{j} = B.\beta = \text{B.B}^{T}.x_{j}$$

Podemos expresar el problema de reducir la dimensionalidad de la
siguiente manera, de acuerdo a lo visto en complemento ortogonal de un
subespacio vectorial:

$${\widetilde{x}}_{j} = \underset{\begin{matrix}
\text{Subespacio} \\
\text{principal} \\
\end{matrix}}{} + \underset{= 0}{}\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (3)$$

Donde el primer término es el subespacio sobre el que vamos a proyectar
los datos y el segundo término (que representa un subespacio
$\mathbb{R}^{D - M}$ dimensional que es el complemento ortogonal del
subespacio principal) para nuestra proyección es cero ya que haremos la
proyección en el subespacio principal.

Definamos una función $\text{J\ }$que represente el error de hacer la
proyección ortogonal de los datos, que como ya dijimos, puede verse como
la sumatoria de los vectores diferencia (representados en el gráfico 5
por las líneas punteadas).

$${J = \frac{1}{N}\sum_{j = 1}^{N}{||x_{j} - {\widetilde{x}}_{j}||}^{2}\backslash n}{\text{donde\ }{||x_{j} - {\widetilde{x}}_{j}||}^{2} = \left( x_{j} - {\widetilde{x}}_{j} \right)^{T}.{(x}_{j} - {\widetilde{x}}_{j})}$$

Lo que buscamos aquí es encontrar un subespacio tal que minimice esta
función de error. Ahora bien, ${\widetilde{x}}_{j}$ va a depender de los
valores que tomen $\beta_{\text{ij}}$ y $b_{i}$ entonces:

$${\frac{\partial J}{\partial\beta_{\text{ij}}} = \frac{\partial J}{\partial{\widetilde{x}}_{j}}.\frac{\partial{\widetilde{x}}_{j}}{\partial\beta_{\text{ij}}} = \underset{\frac{\partial J}{\partial{\widetilde{x}}_{j}}}{}\text{.\ }\underset{\frac{\partial{\widetilde{x}}_{j}}{\partial\beta_{\text{ij}}}}{} = 0\backslash n}{\Rightarrow \ \frac{\partial J}{\partial\beta_{\text{ij}}} = \frac{- 2}{N}\left( x_{j} - \ \ \sum_{h = 1}^{M}{\beta_{\text{hj}}.b_{h}} \right)^{T}\text{.\ }b_{i} = 0}$$

$$\Rightarrow ({x_{j}}^{T}.b_{i}) - \left( \ \sum_{h = 1}^{M}{\beta_{\text{hj}}.b_{h}} \right)^{T}\text{.\ }b_{i} = 0$$

$$\Rightarrow ({x_{j}}^{T}.b_{i}) - \beta_{\text{ij}}.\underset{= 1}{} = 0$$

$$\Rightarrow \beta_{\text{ij}} = \ {x_{j}}^{T}.b_{i}$$

Donde el segundo término sale de derivar (3) respecto a
$\beta_{\text{ij}}$ suponiendo M=1, es decir para la proyección
unidimensional sobre el subespacio generado por el vector $b_{i}$. En el
segundo paso lo que se hace es usar (3) donde ${\widetilde{x}}_{j}$ es
la proyección de $x_{j}$ en el subespacio principal, por eso la
sumatoria solo llega hasta M. En la línea cuatro se simplificó la
sumatoria ya que el producto punto resultante va a ser cero para todos
los términos donde $i \neq h$ por la propiedad vista de bases
ortonormales, y también por esto es que podemos decir que
${b_{i}}^{T}\text{.\ }b_{i}$=1.

Observando la expresión a la que llegamos vemos que coincide con (2) que
era la que habíamos obtenido para proyecciones ortogonales. Entonces, lo
que estamos diciendo aquí es simplemente que las coordenadas
$\beta_{\text{ij}}$ que minimizan el error J son las que nos proyectan
ortogonalmente el vector $x_{j}$ en el eje generado por $b_{i}$.

Visto en términos del gráfico 5 queremos decidir la pendiente de la
recta U. Para simplificar la resolución matemática primero re
expresaremos la función J para dejarla en función de $b_{i}$. Partiendo
de la definición de ${\widetilde{x}}_{j}$ en (3) y reemplazando
$\beta_{\text{ij}}$ de acuerdo a (2) podemos expresar:

$${\widetilde{x}}_{j} = \sum_{i = 1}^{M}{\beta_{\text{ij}}.b_{i}}$$

$${\widetilde{x}}_{j} = \sum_{i = 1}^{M}{{{(x}_{j}}^{T}.b_{i}).b_{i}}$$

Teniendo en cuenta esta última expresión obtenida, y expresando $x_{j}$
de la ecuación (1) como la suma del subespacio principal más el
complemento ortogonal:

$$x_{j} = \sum_{i = 1}^{M}{{{(x}_{j}}^{T}.b_{i}).b_{i}}\  + \ \sum_{i = M + 1}^{D}{{{(x}_{j}}^{T}.b_{i}).b_{i}}$$

Podemos entonces re escribir nuestro vector diferencia como:

$$x_{j} - {\widetilde{x}}_{j} = \ \sum_{i = M + 1}^{D}{{{(x}_{j}}^{T}.b_{i}).b_{i}}$$

Reemplazando esto en nuestra función J usando el hecho de que el
producto punto ${x_{j}}^{T}.b_{i} = {\text{\ b}_{i}}^{T}.x_{j}$ se puede
expresar vectorialmente de cualquiera de esas dos formas:

$$J = \frac{1}{N}\sum_{j = 1}^{N}\left\| \sum_{i = M + 1}^{D}{({\text{\ b}_{i}}^{T}.x_{j}).b_{i}} \right\|^{2}$$

$$\Rightarrow J = \frac{1}{N}\sum_{j = 1}^{N}{\sum_{i = M + 1}^{D}{\left\| {\text{\ b}_{i}}^{T}.x_{j} \right\|^{2}.\underset{= 1}{}}}$$

$$\Rightarrow J = \frac{1}{N}\sum_{j = 1}^{N}{\sum_{i = M + 1}^{D}{({\text{\ b}_{i}}^{T}.x_{j}}).({x_{j}}^{T}.b_{i})}$$

$$\Rightarrow J = \sum_{i = M + 1}^{D}{{\text{\ b}_{i}}^{T}.}\underset{S}{}.b_{i}\ \ \ \ \ \ \ \ \ \ \ \ \ \ (4)$$

Expliquemos brevemente el proceso matemático precedente. En la segunda
línea se saca la sumatoria fuera de la norma ya que va a ser lo mismo
sacar la suma de los vectores y sacarle la norma que sacar la norma de
todos los vectores y sumar dichas normas, también se hizo uso de que la
norma de cada vector base es igual a 1 por bases ortonormales. En la
tercera línea se usó la definición de norma para expresarla como
producto punto. En la cuarta línea se saca “afuera” la sumatoria que
depende de *i* y quedan dentro de la sumatoria con subíndice j solo los
términos que dependen de j.

Analizando la matriz que llamamos S, va a ser la suma de N matrices
$D \times D$. Por lo que cada elemento de la matriz resultante va a ser
la suma de ese elemento de cada una de las N matrices de la sumatoria.
Adicionalmente tenemos que multiplicar cada uno de esos elementos de la
matriz por $\frac{1}{N}$. Podemos ver entonces, que ya que los datos
tienen media igual a cero, esta matriz es la matriz de varianzas y
covarianzas ($D \times D$) de las D variables iniciales. Si los datos
estuvieran estandarizados (adicionalmente estuvieran divididos por su
desviación estándar) la matriz que obtendríamos aquí sería la matriz de
correlaciones.

Supongamos el caso en que tenemos dos variables iniciales y queremos
reducir los datos a una sola componente principal
($\mathbb{R}^{2}\mathbb{\rightarrow R}$). Si planteamos un lagrangiano
para obtener la condición de minimización a partir de la fórmula
obtenida en (4) y teniendo en cuenta la restricción de bases
ortonormales quedará entonces:

$$J = {b_{2}}^{T}\text{.S.}b_{2}\text{\ \ \ \ \ \ \ \ \ \ \ \ s.a.\ }{\text{\ \ \ b}_{2}}^{T}.b_{2} = 1$$

$$L = {b_{2}}^{T}\text{.S.}b_{2} + \lambda(1 - {b_{2}}^{T}.b_{2})$$

$$\frac{\partial L}{\partial\lambda} = 1 - {b_{2}}^{T}.b_{2} = 0$$

$$\frac{\partial L}{\partial b_{2}} = {{2.b}_{2}}^{T}.S - 2.\lambda.{b_{2}}^{T} = 0$$

$$\Rightarrow {b_{2}}^{T}.S = \lambda.{b_{2}}^{T}$$

$$\Rightarrow S.b_{2} = \lambda.b_{2}\ \ \ \ \ \ \ \ \ \ (5)$$

$$J = {b_{2}}^{T}.b_{2}\ .\lambda = \lambda$$

Donde en la cuarta línea se usó la regla de derivada multivariante para
derivar respecto a $b_{2}$. En la sexta línea por ser S una matriz
simétrica puede conmutar de esa forma la multiplicación. En tanto que en
la última línea solo reemplaza el resultado anterior en la función J.

Analicemos la expresión (5), lo que tenemos ahí es un problema de
eigenvalores y eigenvectores. Para recordar, dada una matriz A
($n \times n$), $v$ será eigenvector de A si se cumple que:

$$A.v = \lambda.v$$

Donde $\lambda$ es un escalar y es el eigenvalor de la matriz A
correspondiente a ese eigenvector. Para minimizar entonces J lo que se
debe hacer es buscar los eigenvectores de la matriz S y elegir como
subespacio principal el eigenvector que tenga el $\lambda$. De esta
manera su complemento ortogonal que será el omitido tendrá el
$J = \lambda$ más chico.

Notar que los eigenvectores serán los vectores base, por lo cual por
bases ortonormales serán ortogonales entre sí, y el eigenvector con el
eigenvalor más alto es un vector que apunta en la dirección en la que
los datos tienen más variación (y el eigenvalor va a ser el valor de esa
variación).

Si ampliamos esta lógica para el caso D-dimensional tendremos que:

$$\text{S.}b_{j} = \lambda_{j}.b_{j}\ \ con\ j = M + 1,M + 2,\ldots,D$$

$$J = \ \sum_{j = M + 1}^{D}\lambda_{j}$$

Es decir, minimizamos J eligiendo los D-M eigenvalores más chicos y, por
lo tanto, proyectando los datos en el subespacio formado por los
vectores correspondientes a los M eigenvalores mayores.

Algo a remarcar en este punto, es el hecho de que vamos a tener D
eigenvalores (ya que la dimensión de la matriz S es DxD), de los cuales
elegiremos los M mayores y descartaremos los D-M más chicos.

Analicemos entonces los resultados obtenidos al proyectar los vectores
ortogonalmente sobre el subespacio principal
$\ {x_{j}}^{T}.b_{i} = \ \beta_{\text{ij}}$. Notando esto matricialmente
para el caso de N vectores (recordemos que tenemos N vectores porque
tenemos N observaciones) tendríamos:

$$\underset{N \times D}{} \times \underset{D \times M}{} = \underset{N \times M}{}$$

$$\begin{bmatrix}
x_{11} & \ldots & x_{D1} \\
 \vdots & \ddots & \vdots \\
x_{1N} & \ldots & x_{\text{DN}} \\
\end{bmatrix} \times \begin{bmatrix}
b_{11} & \ldots & b_{M1} \\
 \vdots & \ddots & \vdots \\
b_{1D} & \ldots & b_{\text{MD}} \\
\end{bmatrix} = \begin{bmatrix}
z_{11} & \ldots & z_{M1} \\
 \vdots & \ddots & \vdots \\
z_{1N} & \ldots & z_{\text{MN}} \\
\end{bmatrix}$$

$$z_{11} = \ x_{11}\text{.\ }b_{11} + \ldots + {X_{D1}\text{.b}}_{1D}$$

$$z_{1N} = \ x_{1N}\text{.\ }b_{11} + \ldots + {X_{\text{DN}}\text{.b}}_{1D}$$

$$Z_{1} = \ X.\ b_{1}$$

Cada una de las N observaciones de la nueva variable $z_{1}$ va a ser
una combinación lineal de las D variables iniciales, combinada usando
como coeficientes los valores del eigenvector $b_{1}$. Cada columna de
la matriz Z va a ser una de las nuevas M variables, que no es ni más ni
menos que una de las M coordenadas vistas desde el punto de vista
geométrico de las proyecciones ortogonales.

Mirando los resultados obtenidos para proyecciones ortogonales
M-dimensionales, tenemos que ${\widetilde{x}}_{j} = B.\beta$. La única
aclaración metodológica que corresponde hacer es que aquí estamos viendo
cada uno de los N vectores horizontalmente como una de las filas de la
matriz $\widetilde{X}$, es decir que cada fila podría ser expresada de
acuerdo a la ecuación de proyecciones ortogonales de la siguiente manera
$\underset{1xD}{} = \underset{1xM}{}.\underset{\text{MxD}}{}$.
Expresando esto matricialmente para el caso de N vectores:

$$\underset{N \times M}{} \times \underset{M \times D}{{}^{T}} = \ \underset{N \times D}{}$$

Esto lo que nos está diciendo es que podemos “reconstruir” nuestros
datos originales usando la matriz B como decodificador, la pérdida de
información que vamos a tener va a ser la que perdimos al omitir el
subespacio D-M dimensional complementario al subespacio principal.

Retomemos la expresión $Z_{1} = \ X.\ b_{1}$ y calculemos la esperanza y
la varianza de nuestra nueva variable:

$$E\left( Z_{1} \right) = E\left( \text{X.\ }b_{1} \right)\  = E\left( X \right)\text{.\ }b_{1} = 0\ $$

$$\text{VAR}\left( Z_{1} \right) = \ \frac{1}{N}\text{\ .}{Z_{1}}^{T}.Z_{1} = \ \frac{1}{N}\text{\ .\ }{b_{1}}^{T}\text{.\ }X^{T}\text{\ .\ X\ .\ }b_{1} = \ {b_{1}}^{T}\text{.\ }\underset{S}{}\text{\ .\ }b_{1}$$

Donde en la primera ecuación se usó el hecho de que nuestros datos
iniciales deben estar centrados respecto a su media.

Si observamos la expresión obtenida para la varianza de $Z_{1}$ vemos
que es igual a la obtenida cuando buscábamos cada una de las D-M
dimensiones omitidas, comparando ambos resultados tenemos que:

$${b_{1}}^{T}.S = \lambda.{b_{1}}^{T}\ \ \  \Rightarrow \ \ \ \ \underset{\text{VAR}\left( Z_{1} \right)}{} = \ \lambda.\underset{= 1}{}\ \ \ \ \ \  \Rightarrow \ \ \ \ \ \ \ \ VAR\left( Z_{1} \right) = \ \lambda$$

Interpretando esta última expresión es que podemos decir que buscamos
retener la mayor cantidad de varianza posible, o que vamos a elegir los
M componentes principales con mayor varianza. Esto es porque $Z_{1}$ es
la proyección ortogonal de nuestros datos sobre el eje generado por el
vector $b_{1}$, y este eje es elegido de tal forma de maximizar la
varianza de nuestros datos sobre él. Este enfoque nos permite
compatibilizar la idea de que al retener la mayor cantidad de varianza
posible estamos manteniendo la mayor cantidad de información, y es
porque al mantener la mayor cantidad de varianza estamos minimizando el
error de proyección de nuestros datos en sobre ese subespacio o eje.

Una medida de la variabilidad original de los datos podría ser la suma
de las D varianzas de las variables originales X. Por teoría de
diagonalización (cuya demostración escapa al alcance de este trabajo)
podemos decir que:

$$\text{Traza}\left( S \right) = Traza(P^{'}.V\ .P) = \sum_{j = 1}^{D}{\ \lambda_{j}}$$

Por otra parte, la traza de S (la suma de todos los elementos de su
diagonal principal) va a ser la suma de las D varianzas iniciales, que
es como dijimos nuestra variabilidad total de los datos iniciales.
Adicionalmente si las variables iniciales están estandarizadas esta
varianza total será igual a D (cada una de las D variables iniciales
tiene var=1). Tenemos entonces aquí una relación importante entre la
varianza original y la varianza retenida después del ACP. Por ejemplo,
si retuviéramos todas las componentes principales (M=D) no perderíamos
varianza, lo que puede ser visto como que no habría perdida alguna de
información.

Podríamos decir que la proporción o porcentaje de variabilidad que
explica cada componente va a ser:

$$\%\ de\ Var\ explicado\ por\ Z_{j} = \ \frac{\lambda_{j}}{Traza(S)}$$

Algo que hasta acá hemos dado por sentado es M, es decir la cantidad de
componentes que vamos a retener. Si bien podríamos tener distintos
criterios de decisión uno de los más usados es retener todas las
componentes cuya varianza sea mayor a 1 ($\lambda_{i} > 1$), es decir
que sea mayor a la varianza de cualquiera de las variables iniciales
(todas igual a 1). La interpretación conceptual de esto es que la nueva
componente va a ser “mejor” que cualquiera de las variables iniciales
por si sola ya que va a contener más información.

**CAPÍTULO III**

**ANÁLISIS FACTORIAL**

El análisis factorial nos da un enfoque diferente del problema de
reducción de dimensionalidad. El planteo del que partimos es el
siguiente, dado un conjunto de D variables iniciales existen un conjunto
M de variables subyacentes o factores que no son observables
directamente sino a través de combinaciones de las variables iniciales.
Lo que se busca es representar nuestros datos originales con este
conjunto de factores, teniendo la menor pérdida de información posible,
y al igual que en análisis de componentes principales queremos que los
factores sean obtenidos de tal forma de que no estén correlacionados
entre sí. Esto nos agrega dos elementos de juicio al análisis realizado
hasta acá, uno es el *principio de parsimonia,* queremos que la cantidad
de factores sea lo menor posible (hasta aquí habíamos tomado M como dado
y no elegimos un criterio con el cual decidir la cantidad de componentes
retenidos). El otro elemento es el *principio de interpretabilidad,* es
decir queremos que los factores obtenidos puedan ser interpretables.

Para aclarar esto pongamos un ejemplo, supongamos que se quiere analizar
la importancia que los consumidores dan a 14 variables que se consideran
relevantes para la compra de un automóvil. Estas variables son:
reparaciones baratas (RB), amplia gama de colores (GC), interior
espacioso (IE), bajo consumo de gasolina (BC), manejabilidad (MA),
aspecto moderno (AM), valor de recompra alto (RA), confortable (CO),
motor potente (MP), aspecto elegante (AE), cómodo de conducir (CC),
atractivo de línea (AL), maletero amplio (MA) y fácil de aparcar (FA).
Se observa que las 14 variables pueden caracterizarse por cuatro
dimensiones subyacentes relacionadas respectivamente con el confort
(factor I), con el coste-eficiencia (factor II), con la elegancia
(factor III) y con el manejo fácil (factor IV) y no observables
directamente. Por lo tanto, en vez de considerar las 14 variables,
simplificaremos las cosas, de forma que sólo cuatro factores deban
considerarse para caracterizar la estructura subyacente de los datos. En
el gráfico 6 se puede ver este análisis.

*Gráfico 6*

![](media/image8.png){width="4.608333333333333in"
height="1.9212248468941382in"}

Planteando matemáticamente el modelo factorial nos quedaría:

$${x_{1} = l_{11}.F_{1} + \ldots + l_{1M}.F_{M} + e_{1}\backslash n}{\  \vdots \backslash n}{x_{D} = l_{D1}.F_{1} + \ldots + l_{\text{DM}}.F_{M}\  + e_{D}}$$

Se supone a los factores comunes $F_{1},F_{2},\ldots\ ,\ F_{M}$ como
variables estandarizadas (media cero y varianza unitaria) y que además
no están correlacionadas entre sí. Se supone también que la matriz de
covarianzas de los factores específicos es una matriz diagonal (factores
únicos incorrelacionados entre sí) y tienen media igual a cero
($E\left\lbrack e \right\rbrack = 0\ ,\ \ E\lbrack Fe'\rbrack = 0$)

Dado que las variables X son variables tipificadas, su matriz de
covarianzas es igual a la matriz de correlación poblacional S, matriz
que puede descomponerse de la siguiente forma:

$$\mathbf{S}\  = \ E\left( XX' \right) = E(\left( LF + e \right)\left( LF + e \right)') = LE\left( FF' \right)L' + E\left( ee' \right) + LE\left( fe' \right) + E\left( ef' \right)L' = LIL' + \Omega + L0 + 0L'\mathbf{= LL' + \Omega}$$

Expresando esto matricialmente obtenemos:

$$\underset{\text{DxD}}{} = \underset{\text{DxM}}{} \times \underset{\text{MxD}}{} + \underset{\text{DxD}}{}$$

Si analizamos el resultado para el primer elemento de S (que es la
varianza de $x_{1}$) tenemos que:

$$\text{VAR}\left( x_{1} \right) = 1\  = \ \underset{{h_{1}}^{2}}{+}{\omega_{1}}^{2}$$

Donde ${h_{1}}^{2}$ es el porcentaje de varianza de la variable $x_{1}$
explicado por los factores comunes, y se llama ***comunalidad*** y
${\omega_{1}}^{2}$ es la parte de la varianza de $x_{1}$ que la explica
su factor especifico y se la llama ***especificidad**.*

Ahora bien, cómo relacionamos este análisis de factores con nuestro
análisis de componentes principales. Partamos de la solución expresada
matricialmente a la que llegamos en ACP. Es importante prestar atención
al cambio en la notación, ya que hasta aquí usábamos notación matricial
donde cada elemento era una de las N observaciones, la transformación
que hacemos aquí es pasar a expresarlo vectorialmente, donde cada
variable $Z_{i}$ y $X_{i}$ es un vector Nx1 que contiene todas las
observaciones. Teniendo esto en cuenta reexpresemos la solución
obtenida:

$$\underset{N \times D}{} \times \underset{D \times M}{} = \underset{N \times M}{}$$

$$\begin{bmatrix}
x_{11} & \ldots & x_{D1} \\
 \vdots & \ddots & \vdots \\
x_{1N} & \ldots & x_{\text{DN}} \\
\end{bmatrix} \times \begin{bmatrix}
b_{11} & \ldots & b_{M1} \\
 \vdots & \ddots & \vdots \\
b_{1D} & \ldots & b_{\text{MD}} \\
\end{bmatrix} = \begin{bmatrix}
z_{11} & \ldots & z_{M1} \\
 \vdots & \ddots & \vdots \\
z_{1N} & \ldots & z_{\text{MN}} \\
\end{bmatrix}$$

$${\overset{Nx1}{\overbrace{Z_{1}}} = \ \overset{Nx1}{\overbrace{X_{1}}}\text{.\ }\overset{1x1}{\overbrace{b_{11}}} + \ldots + {X_{D}\text{.b}}_{1D}\backslash n}{\  \vdots}$$

$$Z_{M} = \ X_{1}\text{.\ }b_{M1} + \ldots + {X_{D}\text{.b}}_{\text{MD}}$$

Ahora bien, también habíamos llegado a la conclusión que
$\underset{N \times M}{} \times \underset{M \times D}{{}^{T}} = \ \underset{N \times D}{}$,
que notado vectorialmente quedaría:

$${{\ \widetilde{X}}_{1} = \ Z_{1}\text{.\ }b_{11} + \ldots + {Z_{M}\text{.b}}_{M1}\ \backslash n}{\  \vdots}$$

$${\widetilde{X}}_{D} = \ Z_{1}\text{.\ }b_{1D} + \ldots + {Z_{M}\text{.b}}_{\text{MD}}$$

Recordemos que la diferencia entre ${\ \widetilde{X}}_{i}$ y $X_{i}$ era
la pérdida de información que teníamos por la reducción de
dimensionalidad, supongamos entonces que realizáramos el ACP sin reducir
dimensionalidad (M=D), es decir conservando los D componentes
principales. Este sistema de ecuaciones quedaría de la siguiente forma:

$${X_{1} = \ Z_{1}\text{.\ }b_{11} + \ldots + {Z_{D}\text{.b}}_{D1}\ \backslash n}{\  \vdots}$$

$$X_{D} = \ Z_{1}\text{.\ }b_{1D} + \ldots + {Z_{D}\text{.b}}_{\text{DD}}$$

Pero el análisis factorial exige que los factores estén estandarizados,
entonces estandarizando nuestras variables (teniendo en cuenta los
resultados obtenidos anteriormente
$E\left( Z_{i} \right) = 0\ y\ var\left( Z_{i} \right) = \lambda_{i}\ $)
y llamando $Y_{i}$ a las variables estandarizadas quedaría:

$$Y_{i} = \ \frac{Z_{i}}{\sqrt{\lambda_{i}}}$$

El cual podría ser reexpresado como:

$${\text{\ X}_{1} = \ \underset{\ {\widetilde{X}}_{1}}{} + \underset{e_{1}}{}\ \backslash n}{\ \  \vdots \ }$$

$$X_{D} = \ \underset{{\widetilde{X}}_{D}}{} + \underset{e_{D}}{}$$

Reexpresando cada $Z_{i} = Y_{i}\text{.\ }\sqrt{\lambda_{i}}$ para cada
una de las D ecuaciones tendremos:

$$X_{j} = \ \underset{{\widetilde{X}}_{j}}{} + \ e_{D}$$

Comparando este resultado con el sistema de ecuaciones planteado al
inicio del análisis factorial:

$${x_{1} = l_{11}.F_{1} + \ldots + l_{1M}.F_{M} + e_{1}\backslash n}{\  \vdots \backslash n}{x_{D} = l_{D1}.F_{1} + \ldots + l_{\text{DM}}.F_{M}\  + e_{D}}$$

Donde cada
$l_{\text{ij}} = \ \sqrt{\lambda_{i}}\text{\ .\ }b_{\text{ij}}\text{\ \ \ \ \ \ \ \ con\ \ \ }i = 1,2,..,M\ \ \ \ \ \ \ \ \ j = 1,2,\ldots,D$

Vemos entonces que estandarizando cada uno de los componentes
principales obtenidos del ACP podemos interpretar los resultados del
mismo como un análisis factorial, obteniendo así la ventaja de que los
factores obtenidos puedan tener una interpretación a partir de la
rotación de las componentes que se analizará en la siguiente sección.

Para finalizar esta sección vamos a demostrar cómo cada uno de los
$l_{\text{ij}}$ puede ser interpretado como la correlación entre la
variable inicial “j” y la componente “i”. Notando vectorialmente $X_{j}$
y $\ Z_{i}$ como los vectores con las N observaciones para la variable
“j” y la componente “i” respectivamente podemos calcular su covarianza
como :

$$\text{cov}\left( X_{j},Z_{i} \right) = \ \frac{1}{N}\ {X_{j}}^{'}\text{.\ }\text{\ Z}_{i}$$

Si agregamos un vector $\delta$ de dimensión $1xD$ que tenga 0 en todas
las posiciones y 1 en la posición “j” podríamos notar a ${X_{j}}^{'}$
como:

$${X_{j}}^{'} = \delta\ .\ X^{'} = \ \left\lbrack \begin{matrix}
0 & \ldots & 1 \\
\end{matrix}\text{\ \ \ \ \ }\begin{matrix}
\ldots & 0 \\
\end{matrix} \right\rbrack\text{\ x\ }\begin{bmatrix}
x_{11} & \ldots & x_{1N} \\
 \vdots & \ddots & \vdots \\
x_{D1} & \ldots & x_{\text{DN}} \\
\end{bmatrix}$$

Siendo $\ X^{'}$ la matriz de datos iniciales traspuesta. Adicionalmente
recordando que $Z_{i} = \ X.\ b_{i}$:

$$\text{cov}\left( X_{j},Z_{i} \right) = \ \frac{1}{N}\text{\ δ\ .\ }X^{'}\text{.\ \ X.\ }b_{i} = \ \underset{S}{\ \ }\text{\ .\ \ }\underset{b_{\text{ij}}}{} = \ \lambda_{i}\text{\ .\ }b_{\text{ij}}$$

$$\text{corr}\left( X_{j},Z_{i} \right) = r_{\text{ij}} = \frac{\lambda_{i}\text{\ .\ }b_{\text{ij}}}{\underset{= 1}{\text{\ .\ }\sqrt{\lambda_{i}}\ }}$$

$$r_{\text{ij}} = \sqrt{\lambda_{i}}\text{\ .\ }b_{\text{ij}}$$

Donde se hace uso de la expresión que habíamos obtenido para la matriz
de correlaciones S y de la definición de coeficiente de correlación para
2 variables. De esta forma llegamos a la demostración de que cada uno de
los coeficientes $l_{\text{ij}}$ representa la correlación entre la
variable inicial “j” y la componente “i”.

La matriz formada por todos los coeficientes $l_{\text{ij}}$ es entonces
la **matriz de cargas factoriales** que contiene cada una de las
correlaciones entre los factores retenidos y las variables iniciales.

**CAPÍTULO IV**

**ROTACIÓN DE COMPONENTES**

En el proceso de rotación de componentes lo que se hace gráficamente es
girar la dirección de estos nuevos ejes de referencia que van a definir
el subespacio sobre el que proyectaremos los datos. Existen 2 tipos de
rotaciones, las rotaciones ortogonales que mantienen los ejes
perpendiculares, con lo cual las componentes resultantes seguirán
teniendo la característica de no estar correlacionadas entre sí; y las
rotaciones oblicuas las cuales sacrifican un poco de esa independencia
con el objetivo de obtener una mayor interpretabilidad de las
componentes. A su vez dentro de cada tipo hay varios métodos que se
pueden usar, nosotros en este apartado nos centraremos en la rotación
**Varimax**, que es una rotación ortogonal y es la de uso más extendido.

Lo que buscamos con la rotación de componentes es que cada una de las
variables tenga una correlación máxima con uno de los factores (es
decir, uno) y cero con el resto de los factores o componentes. De tal
forma que esto facilite la tarea de asociar a cada factor como esa
“variable subyacente” representada por un grupo de las variables
iniciales.

Imaginemos que tenemos un conjunto de datos en $\mathbb{R}^{3}$ (es
decir 3 variables iniciales) que podríamos agrupar dentro de un
determinado elipsoide, y que vamos a representar estos datos en un plano
$\mathbb{R}^{2}$ (es decir 2 componentes). La metodología hasta aquí
usada por el ACP es ir eligiendo iterativamente como ejes de cada
componente el eje sobre el que el conjunto de datos tenga mayor
varianza. En este caso como vamos a tener los datos representados en un
plano, tendremos que nuestros nuevos 2 ejes serán los 2 ejes sobre los
que los datos tengan máxima varianza, y que van a coincidir con los 2
ejes más grandes de este elipsoide imaginario. Lo que hacemos al rotar
los componentes es rotar estos 2 ejes (manteniéndolos perpendiculares
entre sí) sin variar el plano que generan. La consecuencia de esto es
que todos los puntos originales se van a seguir representando sobre el
mismo plano, pero respecto a unos ejes diferentes, es decir solo van a
cambiar sus coordenadas. Esto permite que al rotar los componentes **no
perdamos nada del total de varianza explicado** por los componentes en
su conjunto. Metodológicamente lo que hacemos es que la información
perdida por la primera componente (ya no está en la dirección en la que
los datos tienen máxima varianza, por lo tanto, va a disminuir ese nivel
de varianza captada) va a ser recogida por la segunda. Podríamos
generalizar esto diciendo que la información perdida por las primeras K
componentes va a ser recogida por las ultimas M-K componentes.

Matemáticamente la rotación Varimax lo que hace es calcular una variable
que llama **simplicidad** (${C_{i}}^{2}$), que es la varianza de los
cuadrados de las cargas factoriales ($r_{\text{ij}}$) para un
determinado factor (i). Adicionalmente lo que se hace comúnmente es
aplicar lo que se llama la Normalización de Kaiser, donde cada
$r_{\text{ij}}$ se divide por la comunalidad de la variable inicial “j”
(${h_{j}}^{2}$). Esto se hace para evitar que las $x$ con mayor
comunalidad tengan más influencia (más peso) en la solución final. Una
vez calculada esta simplicidad normalizada para cada uno de los factores
lo que se hace es **maximizar la sumatoria de todas estas
simplicidades** ($C^{2}$). Matemáticamente la expresión que se maximiza
es:

$$\text{Max\ }C^{2} = \sum_{i = 1}^{M}{C_{i}}^{2}\text{\ \ \ \ \ \ \ \ \ \ \ con\ \ \ \ }{C_{i}}^{2} = \frac{1}{D}\ \sum_{j = 1}^{D}{\left( \frac{{r_{\text{ij}}}^{2}}{{h_{j}}^{2}} \right)^{2} - \ }\left( \frac{1}{D}\ \sum_{j = 1}^{D}{\frac{{r_{\text{ij}}}^{2}}{{h_{j}}^{2}}\ } \right)^{2}$$

Dos notas interesantes de casos extremos que se pueden hacer sobre este
tema es que**, si no se reduce la cantidad de variables** (M=D) se
pueden rotar los ejes como se quiera sin perder información, solo
cambian las coordenadas (coeficientes). Y que justamente por la forma
iterativa en que se realiza si tuviéramos un solo componente no
podríamos realizar ninguna rotación, como mínimo se deben tener dos.

**CAPÍTULO V**

**CASO PRÁCTICO**

En esta sección desarrollaremos un ejemplo aplicado del ACP en el que a
partir de una serie de variables iniciales (las cuales se detallan en la
tabla 1) que se consideran representativas para la actividad económica
de Mendoza se intentará obtener una cantidad reducida de componentes que
las representen y a partir de los cuales se pueda construir un índice de
actividad económica para la provincia.

La utilidad práctica que se busca con este caso es que, dado que las
variables usadas están disponibles con periodicidad mensual el índice
resultante también lo será, con lo cual podría servir como un estimador
de la actividad económica de la provincia, que hoy solo existe
anualmente.

*Tabla 1*

  Variable                Descripción
  ----------------------- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  AUTOMOTRIZ              Ventas Mensuales de Automotores Cero Km por Segmento. Mendoza. Enero 2010 - Noviembre 2018
  EELECTRESIDENCIAL       Consumo de Energía Eléctrica Residencial en MWh. Años 2004 - 2019
  EELECTCOMERCIAL         Consumo de Energía Eléctrica General/Comercial en MWh. Años 2004 - 2019
  EELECTINDUSTRIAL        Consumo de Energía Eléctrica Grandes Demandas/Industrial en MWh. Años 2004 - 2019
  ENARGAS TOTAL SISTEMA   Total - En miles de m^3^ de 9300 kcal y en porcentaje. Años 2004-2018
  HOTELMDZ                Demanda hotelera por mes y condición de residencia de los viajeros hospedados. Ciudad de Mendoza. Años 2008-2018
  INDUSTRIA               Índice de ventas industriales a valores constantes y variación porcentual.
  INMUEBLES               Total de inmuebles involucrados en operaciones en el Registro Público de la Propiedad, a través de escrituras públicas (corresponden a la 1ª,3ª y 4ª Circunscripción Judicial) y variación porcentual. Mendoza. Enero 2006-Noviembre 2018
  PATVEHI                 Patentamiento de Vehículos (Autos, Motos y Maquinarias Agrícolas)
  SHOPPING                Índice mensual de ventas de mercaderías y servicios en centros de compras a valores corrientes. Año base 2010. Mendoza. Enero 2010-Noviembre 2018
  SUPERMERCADOS           Ventas a precios constantes por grupo de artículos, en pesos de 2004. Mendoza. Años 2010-2018
  VINO                    Despachos de vinos autorizados para ser liberados al consumo. Mendoza. Años 2004 - 2017

Inicialmente usando el software SPSS se realizó un ACP sobre este
conjunto de variables obteniendo los resultados que se muestran a
continuación (tabla 2). Una aclaración es que en todos los casos se
utilizó como metodología para rotar los ejes la rotación Varimax.

*Tabla 2*

Total Variance Explained

  ------------- ----------------------- --------------------------------------- ------------------------------------- --------- ----------------- ---------------- --------- ----------------- ----------------
  > Component   > Initial Eigenvalues   > Extraction Sums of Squared Loadings   > Rotation Sums of Squared Loadings
                > Total                 > % of Variance                         > Cumulative %
  > 1           > 2.816                 > 23.469                                > 23.469
  > 2           > 2.447                 > 20.393                                > 43.862
  > 3           > 1.753                 > 14.608                                > 58.470
  > 4           > 1.500                 > 12.502                                > 70.971
  > 5           > 1.056                 > 8.801                                 > 79.773
  > 6           > .706                  > 5.884                                 > 85.656
  > 7           > .491                  > 4.095                                 > 89.752
  > 8           > .403                  > 3.360                                 > 93.112
  > 9           > .344                  > 2.868                                 > 95.979
  > 10          > .207                  > 1.722                                 > 97.701
  > 11          > .177                  > 1.472                                 > 99.173
  > 12          > .099                  > .827                                  > 100.000
  ------------- ----------------------- --------------------------------------- ------------------------------------- --------- ----------------- ---------------- --------- ----------------- ----------------

  ---------------------------------------------------- ------------- -------------- --------- --------- ---------
  Rotated Component Matrix
  > AUTOS
  > EELECTCOMER
  > EELECTIND
  > EELECTRES
  > GASTOTAL
  > HOTELMDZ
  > INDUSTRIA
  > INMUEBLES
  > PATVEHI
  > SHOPPING
  > SUPERMERCADOS
  > VINO
  Communalities
  > AUTOS
  > EELECTCOMER
  > EELECTIND
  > EELECTRES
  > GASTOTAL
  > HOTELMDZ
  > INDUSTRIA
  > INMUEBLES
  > PATVEHI
  > SHOPPING
  > SUPERMERCADOS
  > VINO
  > Extraction Method: Principal Component Analysis.
  ---------------------------------------------------- ------------- -------------- --------- --------- ---------

Dados que estos resultados nos dejan 5 componentes principales, se
complicaría obtener un índice único. El proceso que se fue haciendo
implicó sacar algunas variables en base al peso que tuvieran en los
componentes principales, descartando las de menor peso relativo. También
en base al número de observaciones disponibles para esa variable, ya que
no de todas se disponía de la misma cantidad de observaciones,
priorizando mantener aquellas que tuvieran más. El resultado final de
este proceso de iteración es el que se muestra en la tabla 3.

*Tabla 3*

Total Variance Explained

  ------------- ----------------------- --------------------------------------- ------------------------------------- --------- ----------------- ---------------- --------- ----------------- ----------------
  > Component   > Initial Eigenvalues   > Extraction Sums of Squared Loadings   > Rotation Sums of Squared Loadings
                > Total                 > % of Variance                         > Cumulative %
  > 1           > 2.403                 > 40.054                                > 40.054
  > 2           > 1.673                 > 27.879                                > 67.933
  > 3           > .881                  > 14.689                                > 82.622
  > 4           > .564                  > 9.398                                 > 92.020
  > 5           > .368                  > 6.132                                 > 98.152
  > 6           > .111                  > 1.848                                 > 100.000
  ------------- ----------------------- --------------------------------------- ------------------------------------- --------- ----------------- ---------------- --------- ----------------- ----------------

  -------------------------- ------------- --------------
  Rotated Component Matrix
  > EELECTCOMER
  > EELECTRES
  > GASTOTAL
  > INMUEBLES
  > VINO
  > HOTELMDZ
  Communalities
  > EELECTCOMER
  > EELECTRES
  > GASTOTAL
  > INMUEBLES
  > VINO
  > HOTELMDZ
  -------------------------- ------------- --------------

De acuerdo a estos resultados finales usamos 6 variables representativas
del nivel de actividad mensual y obtenemos 2 factores principales. De
las 6 variables EELECTCOMER, EELECTRES y HOTELMDZ están muy
representadas en el factor 1, en tanto que GASTOTAL, INMUEBLES y VINO lo
están en el factor 2.

Dados estos 2 factores obtenidos, cuyos histogramas están representados
en el gráfico 7, los combinamos algebraicamente sumándolos, para así
obtener el índice de actividad que estábamos buscando.

Ahora bien, para probar la validez de nuestro índice vamos a compararlo
respecto al EMAE en el mismo periodo de tiempo, para que esta
comparación tenga sentido lo que se hizo fue estandarizar el EMAE como
se muestra en la Tabla 4 y Gráfico 8 respectivamente

*Gráfico 7*

![](media/image9.png){width="3.5027777777777778in"
height="2.738888888888889in"}![](media/image10.png){width="4.676388888888889in"
height="2.7527777777777778in"}

*Tabla 4*

  ------------------------------ ------- ----------- ----------- ------------ ------------------
  > **Descriptive Statistics**
  > emae2
  > Valid N (listwise)
  ------------------------------ ------- ----------- ----------- ------------ ------------------

*Gráfico 8*

![](media/image11.png){width="6.6in" height="3.885255905511811in"}

Luego de estandarizado el EMAE lo contraponemos gráficamente con
respecto a ambos factores (gráfico 9) y con respecto a nuestro índice
(gráfico 10).

*Gráfico 9*

[\[CHART\]]{.chart}

*Gráfico 10*

[\[CHART\]]{.chart}

Como podemos observar en el gráfico 11, sí hay una equivalencia en las
variaciones de los niveles de actividad a nivel nacional con respecto a
lo que muestra el índice construido para Mendoza. Adicionalmente se
puede observar un aparente rezago del índice respecto del EMAE, para
corroborar esto se intentará ver si este rezago existe también entre el
PBG de la provincia y el PBI nacional que serían los equivalentes de
estas series, pero con periodicidad anual. Esto último se muestra en los
gráficos 11 y 12 respectivamente.

*Gráfico 11*

[\[CHART\]\[CHART\]]{.chart}

Analizando los resultados obtenidos de nuestro caso de estudio, podemos
concluir que a través del método de ACP fuimos capaces de crear un
índice que representa la actividad económica de la provincia. Al
comparar este índice respecto al EMAE hay una equivalencia con algunos
períodos de atraso, lo cual indicaría que las variaciones que se
producen a nivel nacional se transmiten luego a la provincia. Al momento
de intentar contrastar esta hipótesis haciendo uso del PBG y PBI, la
misma pareciera no cumplirse (teniendo en cuenta que la periodicidad es
distinta en ambos casos ya que justamente el índice se crea porque no
hay ninguna estadística con esa periodicidad para la provincia) sino que
son bastante simultaneas las variaciones. Quedaría como posible caso de
estudio profundizar en el por qué de este rezago al hacer el análisis
mensual.

CONCLUSIONES

Para sacar las conclusiones usaremos ejemplos de
$\mathbb{R}^{2}\mathbb{\rightarrow R}$ para que sea más simple
analizarlo conceptualmente más allá de la matemática. Supongamos que
tenemos una serie de datos con una representación gráfica como la del
gráfico 12.

*Gráfico 12*

![D:\\Asset 6.png](media/image12.png){width="3.6402012248468942in"
height="3.460870516185477in"}

Calculamos los eigenvectores y el resultado obtenido nos dice que la
proyección que minimiza la función J es la dada por la recta $y = 2x$. O
sea que nuestros errores van a ser las distancias de cada punto a esa
recta. Entre más correlacionadas estuvieran estas dos variables esta
nube de puntos más se parecería a una recta y por ende menor sería el
error. Esto nos da el primer caso en el que puede ser útil el ACP,
cuando las variables iniciales están altamente correlacionadas podemos
reducir su dimensionalidad sin tener una gran pérdida de información.

Otro caso sería que, por ejemplo, la nube de puntos estuviera alrededor
del eje X, es decir la recta tendría pendiente igual a cero. En este
caso lo que está pasando es que la variable Y tiene muy poca varianza.
El resultado del ACP para este caso va a ser que el subespacio principal
va a ser el eje X, cuyo vector base es $b = \begin{bmatrix}
1 \\
0 \\
\end{bmatrix}$. Esto quiere decir que el coeficiente que le estamos
asignando a los datos de la variable Y es 0, estamos eliminando esta
variable. Este es el segundo caso, cuando haya variables con varianza
muy bajas van a tener coeficientes muy cercanos a cero y
$J\widetilde{=}\ \text{var}_{y}$ (en nuestro caso extremo J es
exactamente igual a la varianza de Y, la perdida que tenemos de
información es la varianza de la variable que estamos eliminando). O sea
que viendo la matriz S a priori nos podríamos dar cuenta que tanta
pérdida de información tendríamos por reducir su dimensionalidad,
valores altos de correlaciones o valores bajos de varianzas indicarían
casos donde se podría aplicar ACP con muy poca pérdida de información.

Una visión alternativa para el uso del ACP sería como método de
compresión de datos. Siguiendo este ejemplo
$\mathbb{R}^{2}\mathbb{\rightarrow R}$ inicialmente tendremos dos series
de N datos, y el ACP nos deja una sola serie de N coordenadas. Si a esa
serie de coordenadas la multiplicamos por el vector base nos da los N
puntos de nuevo bidimensionales, pero ahora todos sobre la recta
$y = 2x$. Visto para grandes cantidades de datos esto nos permitiría
almacenarlos solo con las series de coordenadas que tienen una
dimensionalidad menor y la matriz de coordenadas haría las veces de
decodificador para volver a obtener la información inicial.

Una última conclusión y no poco importante es que la componente
principal obtenida (la “nueva variable”) va a ser una combinación lineal
de las variables iniciales. El vector base que genera la recta $y = 2x$
es $b = \begin{bmatrix}
0.45 \\
0.9 \\
\end{bmatrix}$ (aproximadamente). Entonces el valor que va a tener la
coordenada, que puede ser vista como una nueva variable Z va a ser los
datos de cada par (X,Y) combinados linealmente con esos coeficientes.
Esto recordemos sale del resultado de la coordenada para proyecciones
ortogonales ($B_{\text{ij}} = {x_{j}}^{T}\text{\ .}b_{i}$)

BIBLIOGRAFÍA CONSULTADA
