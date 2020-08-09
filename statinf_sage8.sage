#!/usr/bin/env sage -python

from sage.all import *

def mostrar_decimal(x, desc='', fmt='{:4.3f}'):
    text = desc + fmt
    print(text.format(float(x)))

def mostrar_numero(x, digitos=3):
    return(x.n(digits=digitos))

class Estadin(object):
    """
    DESCRIPCIÓN
    Clase con funciones para cálculos de inferencia
    estadística.
    """
    def pnorm(self, z):
        """
        Obtener la probabilidad acumulativa
        normal o gaussiana para el cuantil z.
        """
        N = RealDistribution('gaussian', 1)
        p = N.cum_distribution_function(z)
        return(p)

    def qnorm(self, p):
        """
        Obtener el cuantil o distribución acumulativa
        inversa normal o gaussiana para la probabilidad p.
        """
        N = RealDistribution('gaussian', 1)
        z = N.cum_distribution_function_inv(p)
        return(z)

    def pt(self, z, df):
        """
        Obtener la probabilidad acumulativa
        t-Student con df grados de libertad.
        """
        T = RealDistribution('t', df)
        p = T.cum_distribution_function(z)
        return(p)

    def qt(self, p, df):
        """
        Obtener el cuantil t-Student con df grados
        de libertad para la probabilidad p.
        """
        T = RealDistribution('t', df)
        z = T.cum_distribution_function_inv(p)
        return(z)

    def pchisq(self, z, df):
        """
        Obtener la probabilidad acumulativa
        Chi-squared con df grados de libertad.
        """
        C = RealDistribution('chisquared', df)
        p = C.cum_distribution_function(z)
        return(p)

    def qchisq(self, p, df):
        """
        Obtener el cuantil Chi-squared con df grados
        de libertad para la probabilidad p.
        """
        C = RealDistribution('chisquared', df)
        z = C.cum_distribution_function_inv(p)
        return(z)

    def obtener_radio_ic_normal(self, confianza):
        """
        Obtener el radio normal para la confianza dada.
        """
        N = RealDistribution('gaussian', 1)
        P = (1+confianza)/2
        r = N.cum_distribution_function_inv(P)
        return(r)

    def obtener_amplitud_muestra_proporcion(
        self, radio_ic, confianza):
        """
        Obtener el amplitud de una muestra para estimar
        una proporción para la confianza dada.

        DESCRIPCIÓN: El amplitud de la muestra es
        directamente proporcional de la desviación
        estándar de la distribución de éxitos de la
        población, que depende de la proporción p de
        éxitos y es la raíz cuadrada de p(1-p).
        Si bien desconocemos p, como es menor o igual
        a 1, p(1-p) es menor o igual a 1/4.
        Entonces el extremo superior de los amplituds
        de muestras es (1/2) (z/radio_ic)^2. Donde z es el
        quantil normal de la confianza.

        ENTRADA:
        radio_ic: radio_ic dada
        confianza: confianza dada
        SALIDA:
        amplitud de la muestra
        TEST:
        >>> cest = Estadin()
        >>> cest.obtener_amplitud_muestra(0.01, 0.95)
        9604
        """
        r = self.obtener_radio_ic_normal(confianza)
        n = 0.25*(r/radio_ic)^2
        return(n.ceil())

    def obtener_ic_proporcion(self, tam, confianza, prop,
                               fmt=':3.2f'):
        """
        Obtener el intervalo de confianza para una
        proporción.

        ENTRADA:
        tam: amplitud de la muestra,
        confianza: confianza deseada,
        prop: proporción de éxitos de la muestra.

        SALIDA:
        (prop-c, prop+c): donde c es el radio
        de confianza, obtenido a partir de los datos
        de entrada.

        DESCRIPCIÓN:
        Asumiendo que prop es una variable aleatoria
        que corresponde a dividir por tam una suma de
        variables aleatorias Bernoulli(p), y por TLC
        tendrá una distribución aproximadamente normal
        de promedio p y desviación estándar
        s = sqrt(prop*(1-prop)/tam).
        Normalizando: z = (prop - p)/s,
        tendrá distribución normal de promedio 0 y
        desviación estándar 1.

        TEST:

        >>> ces = Estadin()
        >>> ces.obtener_ic_proporcion(tam = 100,
                                       confianza = 0.95,
                                       prop = 35/100,
                                       fmt=':3.2f')
        '(0.26, 0.44)'

        """
        r = self.obtener_radio_ic_normal(confianza)
        var = prop*(1-prop)/tam
        s = sqrt(var)
        radio_ic = r*s
        ic = (prop-radio_ic, prop+radio_ic)
        res = '({'+fmt+'}, {'+fmt+'})'
        ica = float(ic[0])
        icb = float(ic[1])
        intervalo = res.format(ica, icb)
        return(intervalo)

    def obtener_ic_promedio(self, tam, confianza, prom,
                            s, fmt=':3.0f'):
        """
        Obtener el intervalo de confianza para un
        promedio.

        ENTRADA:
        tam: amplitud de la muestra,
        confianza: confianza deseada,
        prom: promedio de valores de la muestra,
        s: desviación estándar:
            s = sqrt(1/(n-1) sum_1^n (X_i - E(X))^2)

        SALIDA:
        (prom-c, prom+c): donde c es el radio
        de confianza, obtenido a partir de los datos
        de entrada.

        DESCRIPCIÓN:
        Asumiendo que prom es una variable aleatoria
        que corresponde a dividir por tam una suma de
        variables aleatorias, y por TLC
        tendrá una distribución aproximadamente normal
        de promedio p y varianza s^2/n.

        Normalizando, la variable aleatoria
        z = (prop - p)/s, tendrá distribución normal
        de promedio 0 y desviación estándar 1.

        TEST:

        >>> ces = Estadin()
        >>> ces.obtener_ic_promedio(tam = 500,
                                    confianza = 0.90,
                                    prom = 562,
                                    s = 112,
                                    fmt=':3.0f')
        '(554, 570)'

        """
        r = self.obtener_radio_ic_normal(confianza)
        # var = prom*(1-prom)/tam
        radio_ic = r * s / sqrt(tam)
        #print('radio_ic = '.format(float(radio_ic)))
        ic = (prom-radio_ic, prom+radio_ic)
        res = '({'+fmt+'}, {'+fmt+'})'
        ica = float(ic[0])
        icb = float(ic[1])
        intervalo = res.format(ica, icb)
        return(intervalo)

    def obtener_ic_dif_proporciones(self, confianza, tam1,
                                    prop1, tam2, prop2,
                                    fmt=':3.2f'):
        """
        Obtener el intervalo de confianza para una
        diferencia de proporciones.

        ENTRADA:
        confianza: confianza deseada,
        tam1: amplitud de la muestra1,
        prop1: proporción de éxitos de la muestra1.
        tam2: amplitud de la muestra2,
        prop2: proporción de éxitos de la muestra2.

        SALIDA:
        (prop1-prop2-c, prop1-prop2+c): donde c es el radio
        de confianza, obtenido a partir de los datos de
        entrada.

        DESCRIPCIÓN:

        TEST:

        >>> ces = Estadin()
        >>> ces.obtener_ic_dif_proporciones(confianza = 0.9,
                                       tam1 = 100,
                                       prop1 = 35/100,
                                       tam2 = 200,
                                       prop2 = 50/200,
                                       fmt=':3.2f')
        '(0.01, 0.19)'
        """
        r = self.obtener_radio_ic_normal(confianza)
        var1 = prop1*(1-prop1)/tam1
        var2 = prop2*(1-prop2)/tam2
        var = var1+var2
        s = sqrt(var)
        radio_ic = r*s
        ic = (prop1-prop2-radio_ic, prop1-prop2+radio_ic)
        res = '({'+fmt+'}, {'+fmt+'})'
        ica = float(ic[0])
        icb = float(ic[1])
        intervalo = res.format(ica, icb)
        return(intervalo)

    def obtener_ic_dif_promedios(self, confianza, tam1,
                                 prom1, s1, tam2, prom2, s2,
                                 fmt=':3.2f'):
        """
        Obtener el intervalo de confianza para una
        diferencia de promedios.

        ENTRADA:
        confianza: confianza deseada,
        tam1: amplitud de la muestra1,
        prom1: promedio de la muestra1,
        s1: desviación estándar de la muestra1,
        tam2: amplitud de la muestra2,
        s2: desviación estándar de la muestra2,
        prop2: proporción de éxitos de la muestra2.

        SALIDA:
        (prom1-prom2-c, prom1-prom2+c): donde c es el radio
        de confianza, obtenido a partir de los datos
        de entrada.

        DESCRIPCIÓN:

        TEST:

        >>> ces = Estadin()
        >>> ces.obtener_ic_dif_promedios(confianza = 0.95,
                                         tam1 = 500,
                                         prom1 = 562,
                                         s1 = 112,
                                         tam2 = 300,
                                         prom2 = 551,
                                         s2 = 121,
                                         fmt=':3.2f')
        '(-6, 28)'
        """
        r = self.obtener_radio_ic_normal(confianza)
        var1 = s1^2/tam1
        var2 = s2^2/tam2
        var = var1+var2
        s = sqrt(var)
        radio_ic = r*s
        ic = (prom1-prom2-radio_ic, prom1-prom2+radio_ic)
        res = '({'+fmt+'}, {'+fmt+'})'
        ica = float(ic[0])
        icb = float(ic[1])
        intervalo = res.format(ica, icb)
        return(intervalo)


    def obtener_puntos_densidad(self, tipod=('gaussian', 1),
                                r=2, zi=-3, zf=4):
        """
        ENTRADA:
        tipod: (nombre, parámetro) de la distribución
        r: 2^r cantidad de puntos por unidad
        zi: extremo inferior del intervalo,
        zf: extremo superior del intervalo
        DESCRIPCIÓN:
        n: total de puntos
        x_i: abscisas de los puntos para cada i,
            en [0,n*2^r]
        """
        test_orden = zi < zf
        if not test_orden:
            print ("INTERVALOS MAL DEFINIDOS")
            return None

        D = RealDistribution(tipod[0],tipod[1])
        n = (zf-zi).ceil()*2^r
        if n <= 0:
            return [(zi,0),(zf,0)]
        x_i= lambda i: zi + i*(zf-zi)/(n*2^(r))
        dpoints = [
            (x_i(i), D.distribution_function(x_i(i)))
            for i in range(0,1+n*2^(r)) if x_i(i) <= zf]
        return([(zi,0)] + dpoints + [(zf,0)])

    def dibujar_densidad(self, tipod=('gaussian',1), r=4,
                         zl=-4, zi=None, zf=None, zr=4,
                         alpha=0.5, rgbcolor='blue',
                         edgecolor='blue', thickness=1,
                         legend_label=None,
                         legend_color=None, aspect_ratio=8,
                         fill=True, complemento=False):
        """
        ENTRADA:
        tipod: (nombre, parámetro) de la distribución
        r: 2^r puntos por unidad
        zl: extremo inferior del intervalo,
        zi,zf: subintervalo,
        zr: extremo superior del intervalo,

        SALIDA:
        Objeto gráfico con el dibujo.

        DESCRIPCIÓN:
        Dibujar la densidad de la distribución
        y de acuerdo al subintervalo, mostrar áreas dentro
        o fuera (complemento = True).
        """
        if zi == None:
            zi=zl
        if zf == None:
            zf = zr

        test_orden = zl <=zi and zi < zf and zf <= zr
        if not test_orden:
            print ("INTERVALOS MAL DEFINIDOS")
            return None

        D = RealDistribution(tipod[0],tipod[1])
        Ptot = D.plot(zl,zr, alpha=alpha, rgbcolor=rgbcolor,
                      thickness=thickness,
                      legend_label=legend_label,
                      legend_color=legend_color,
                      aspect_ratio=aspect_ratio)

        poligono = lambda zi, zf: polygon(
            self.obtener_puntos_densidad(tipod, r, zi, zf),
            alpha=alpha, rgbcolor=rgbcolor,
            edgecolor=edgecolor, thickness=thickness,
            aspect_ratio=aspect_ratio, fill=fill)
        Plot = Ptot
        if not complemento:
            Plot += poligono(zi,zf)
        else:
            if zl < zi:
                Plot += poligono(zl,zi)
            if zf < zr:
                Plot += poligono(zf,zr)
        return(Plot)
