
/* =========================================
   EULER EN EL ESPACIO
   Simulador de lanzamiento de un cohete
   Métodos: Euler, Heun, RK4 y solución analítica
========================================= */


// =========================================
// 1. CONFIGURACIÓN INICIAL
// =========================================

const valoresIniciales = {
    masa: 1000,
    empuje: 15000,
    gravedad: 9.81,
    resistencia: 100,
    velocidadInicial: 0,
    alturaInicial: 0,
    tiempoTotal: 20,
    paso: 2
};

let resultadosActuales = [];
let parametrosActuales = null;


// Formato numérico utilizado en la página.
const formatoNumero = new Intl.NumberFormat("es-BO", {
    minimumFractionDigits: 4,
    maximumFractionDigits: 4,
    useGrouping: false
});


function formato(valor) {
    return formatoNumero.format(valor);
}


// =========================================
// 2. OBTENER LOS DATOS DEL FORMULARIO
// =========================================

function leerParametros() {

    return {
        masa: document.getElementById("mass").valueAsNumber,

        empuje: document.getElementById("thrust").valueAsNumber,

        gravedad: document.getElementById("gravity").valueAsNumber,

        resistencia: document.getElementById("resistance").valueAsNumber,

        velocidadInicial: document.getElementById("initial-velocity").valueAsNumber,

        alturaInicial: document.getElementById("initial-height").valueAsNumber,

        tiempoTotal: document.getElementById("total-time").valueAsNumber,

        paso: document.getElementById("step-size").valueAsNumber
    };
}


// =========================================
// 3. VALIDAR LOS PARÁMETROS
// =========================================

function validarParametros(p) {

    if (Object.values(p).some(valor => !Number.isFinite(valor))) {
        return "Completa todos los campos con valores numéricos válidos.";
    }

    if (p.masa <= 0) {
        return "La masa del cohete debe ser mayor que cero.";
    }

    if (p.tiempoTotal <= 0) {
        return "El tiempo de simulación debe ser mayor que cero.";
    }

    if (p.paso <= 0) {
        return "El paso de integración debe ser mayor que cero.";
    }

    if (p.empuje < 0 || p.gravedad < 0 || p.resistencia < 0) {
        return "El empuje, la gravedad y la resistencia no pueden ser negativos.";
    }

    const cantidadPasos = Math.ceil(p.tiempoTotal / p.paso);

    if (cantidadPasos > 5000) {
        return "El paso seleccionado genera demasiadas iteraciones. Utiliza un paso mayor.";
    }

    return "";
}


// =========================================
// 4. ECUACIÓN DIFERENCIAL DE VELOCIDAD
// =========================================

function aceleracion(velocidad, p) {

    return (p.empuje / p.masa)
        - p.gravedad
        - (p.resistencia / p.masa) * velocidad;
}


// =========================================
// 5. SOLUCIÓN ANALÍTICA
// =========================================

function solucionAnalitica(t, p) {

    const a = (p.empuje / p.masa) - p.gravedad;

    const lambda = p.resistencia / p.masa;

    let velocidad;
    let altura;


    // Caso especial: resistencia igual a cero.
    if (Math.abs(lambda) < 1e-12) {

        velocidad = p.velocidadInicial + a * t;

        altura = p.alturaInicial
            + p.velocidadInicial * t
            + 0.5 * a * t * t;

    } else {

        const velocidadEquilibrio = a / lambda;

        const exponencial = Math.exp(-lambda * t);

        velocidad = velocidadEquilibrio
            + (p.velocidadInicial - velocidadEquilibrio) * exponencial;


        altura = p.alturaInicial
            + velocidadEquilibrio * t
            + (p.velocidadInicial - velocidadEquilibrio)
            * (1 - exponencial) / lambda;

    }

    return {
        velocidad: velocidad,
        altura: altura
    };
}


// =========================================
// 6. MÉTODO DE EULER
// =========================================

function calcularEuler(velocidad, altura, paso, p) {

    // Cálculo de la nueva velocidad.
    const nuevaVelocidad = velocidad
        + paso * aceleracion(velocidad, p);

    // Cálculo de la nueva altura.
    const nuevaAltura = altura
        + paso * velocidad;

    return {
        velocidad: nuevaVelocidad,
        altura: nuevaAltura
    };
}


// =========================================
// 7. MÉTODO DE HEUN
// =========================================

function calcularHeun(velocidad, altura, paso, p) {

    // Primera pendiente.
    const k1v = aceleracion(velocidad, p);

    const k1y = velocidad;


    // Predicción de los valores.
    const velocidadPredicha = velocidad + paso * k1v;


    // Segunda pendiente.
    const k2v = aceleracion(velocidadPredicha, p);

    const k2y = velocidadPredicha;


    // Promedio de las pendientes.
    const nuevaVelocidad = velocidad
        + (paso / 2) * (k1v + k2v);

    const nuevaAltura = altura
        + (paso / 2) * (k1y + k2y);


    return {
        velocidad: nuevaVelocidad,
        altura: nuevaAltura
    };
}


// =========================================
// 8. MÉTODO RUNGE-KUTTA DE CUARTO ORDEN
// =========================================

function calcularRK4(velocidad, altura, paso, p) {

    // Primera evaluación.
    const k1v = aceleracion(velocidad, p);

    const k1y = velocidad;


    // Segunda evaluación.
    const v2 = velocidad + (paso / 2) * k1v;

    const k2v = aceleracion(v2, p);

    const k2y = v2;


    // Tercera evaluación.
    const v3 = velocidad + (paso / 2) * k2v;

    const k3v = aceleracion(v3, p);

    const k3y = v3;


    // Cuarta evaluación.
    const v4 = velocidad + paso * k3v;

    const k4v = aceleracion(v4, p);

    const k4y = v4;


    // Resultado final de velocidad.
    const nuevaVelocidad = velocidad
        + (paso / 6) * (
            k1v + 2 * k2v + 2 * k3v + k4v
        );


    // Resultado final de altura.
    const nuevaAltura = altura
        + (paso / 6) * (
            k1y + 2 * k2y + 2 * k3y + k4y
        );


    return {
        velocidad: nuevaVelocidad,
        altura: nuevaAltura
    };
}


// =========================================
// 9. GUARDAR UNA ITERACIÓN
// =========================================

function crearFila(t, euler, heun, rk4, p) {

    const exacta = solucionAnalitica(t, p);

    return {

        tiempo: t,

        eulerV: euler.velocidad,
        eulerY: euler.altura,

        heunV: heun.velocidad,
        heunY: heun.altura,

        rk4V: rk4.velocidad,
        rk4Y: rk4.altura,

        exactV: exacta.velocidad,
        exactY: exacta.altura

    };
}


// =========================================
// 10. EJECUTAR LOS CUATRO MÉTODOS
// =========================================

function ejecutarSimulacion(p) {

    let tiempo = 0;


    // Estados iniciales de cada método.
    let euler = {
        velocidad: p.velocidadInicial,
        altura: p.alturaInicial
    };

    let heun = {
        velocidad: p.velocidadInicial,
        altura: p.alturaInicial
    };

    let rk4 = {
        velocidad: p.velocidadInicial,
        altura: p.alturaInicial
    };


    const filas = [];


    // Guardar las condiciones iniciales.
    filas.push(crearFila(tiempo, euler, heun, rk4, p));


    // Repetir los cálculos hasta llegar al tiempo final.
    while (tiempo < p.tiempoTotal - 1e-12) {

        // Evitar superar el tiempo final.
        const pasoActual = Math.min(
            p.paso,
            p.tiempoTotal - tiempo
        );


        // Calcular Euler.
        euler = calcularEuler(
            euler.velocidad,
            euler.altura,
            pasoActual,
            p
        );


        // Calcular Heun.
        heun = calcularHeun(
            heun.velocidad,
            heun.altura,
            pasoActual,
            p
        );


        // Calcular RK4.
        rk4 = calcularRK4(
            rk4.velocidad,
            rk4.altura,
            pasoActual,
            p
        );


        // Actualizar el tiempo.
        tiempo += pasoActual;

        if (Math.abs(tiempo - p.tiempoTotal) < 1e-10) {
            tiempo = p.tiempoTotal;
        }


        // Guardar todos los resultados.
        filas.push(
            crearFila(tiempo, euler, heun, rk4, p)
        );

    }

    return filas;
}


// =========================================
// 11. MOSTRAR LAS TARJETAS DE RESULTADOS
// =========================================

function mostrarResultados(filas) {

    const ultimo = filas[filas.length - 1];


    const errorVelocidad = Math.abs(
        ultimo.exactV - ultimo.eulerV
    );

    const errorAltura = Math.abs(
        ultimo.exactY - ultimo.eulerY
    );


    // Resultados principales de Euler.
    document.getElementById("velocity-result").textContent =
        formato(ultimo.eulerV);

    document.getElementById("height-result").textContent =
        formato(ultimo.eulerY);

    document.getElementById("velocity-error").textContent =
        formato(errorVelocidad);

    document.getElementById("height-error").textContent =
        formato(errorAltura);


    mostrarPrecision(errorVelocidad, errorAltura);

    mostrarComparacionMetodos(ultimo);

}


// =========================================
// 12. COMPROBAR LA PRECISIÓN
// =========================================

function mostrarPrecision(errorVelocidad, errorAltura) {

    const banner = document.getElementById("precision-banner");

    const titulo = document.getElementById("precision-title");

    const descripcion = document.getElementById("precision-description");

    const icono = banner.querySelector(".precision-icon");


    const cumpleVelocidad = errorVelocidad < 0.5;

    const cumpleAltura = errorAltura < 5;


    if (cumpleVelocidad && cumpleAltura) {

        banner.className = "precision-banner precision-ok";

        icono.textContent = "✓";

        titulo.textContent = "Precisión alcanzada";

        descripcion.textContent =
            "Los errores de velocidad y altura cumplen los límites establecidos.";

    } else {

        banner.className = "precision-banner precision-warning";

        icono.textContent = "!";

        titulo.textContent = "Precisión insuficiente";

        descripcion.textContent =
            "Reduce el tamaño del paso para mejorar los resultados de Euler.";

    }

}


// =========================================
// 13. COMPARAR LOS MÉTODOS
// =========================================

function mostrarComparacionMetodos(ultimo) {

    const contenedor = document.getElementById("method-cards");


    const metodos = [

        {
            nombre: "Euler explícito",
            clase: "euler",
            velocidad: ultimo.eulerV,
            altura: ultimo.eulerY
        },

        {
            nombre: "Heun · Orden 2",
            clase: "heun",
            velocidad: ultimo.heunV,
            altura: ultimo.heunY
        },

        {
            nombre: "Runge-Kutta · RK4",
            clase: "rk4",
            velocidad: ultimo.rk4V,
            altura: ultimo.rk4Y
        },

        {
            nombre: "Solución analítica",
            clase: "exact",
            velocidad: ultimo.exactV,
            altura: ultimo.exactY
        }

    ];


    contenedor.innerHTML = "";


    metodos.forEach(metodo => {

        const errorV = Math.abs(
            ultimo.exactV - metodo.velocidad
        );

        const errorY = Math.abs(
            ultimo.exactY - metodo.altura
        );


        const tarjeta = document.createElement("article");

        tarjeta.className = "method-card " + metodo.clase;


        const nombre = document.createElement("h3");
        nombre.textContent = metodo.nombre;

        tarjeta.appendChild(nombre);


        const datos = [
            ["Velocidad", formato(metodo.velocidad) + " m/s"],
            ["Altura", formato(metodo.altura) + " m"]
        ];


        datos.forEach(dato => {

            const fila = document.createElement("div");

            fila.className = "method-value";

            const etiqueta = document.createElement("span");
            etiqueta.textContent = dato[0];

            const valor = document.createElement("strong");
            valor.textContent = dato[1];

            fila.append(etiqueta, valor);

            tarjeta.appendChild(fila);

        });


        const errores = document.createElement("div");

        errores.className = "method-error";

        if (metodo.clase === "exact") {

            errores.textContent = "Referencia exacta · Error 0";

        } else {

            errores.textContent =
                "Errores: V = " + formato(errorV)
                + " m/s · Y = " + formato(errorY) + " m";

        }


        tarjeta.appendChild(errores);

        contenedor.appendChild(tarjeta);

    });

}


// =========================================
// 14. MOSTRAR LA TABLA DE ITERACIONES
// =========================================

function mostrarTabla(filas) {

    const tbody = document.getElementById("results-table");

    tbody.innerHTML = "";


    filas.forEach(fila => {

        const tr = document.createElement("tr");


        const valores = [

            formato(fila.tiempo),

            formato(fila.eulerV),
            formato(fila.eulerY),

            formato(fila.heunV),
            formato(fila.heunY),

            formato(fila.rk4V),
            formato(fila.rk4Y),

            formato(fila.exactV),
            formato(fila.exactY)

        ];


        valores.forEach(valor => {

            const td = document.createElement("td");

            td.textContent = valor;

            tr.appendChild(td);

        });


        tbody.appendChild(tr);

    });

}


// =========================================
// 15. DIBUJAR LAS GRÁFICAS
// =========================================

function dibujarGrafica(idCanvas, propiedad, filas) {

    const canvas = document.getElementById(idCanvas);

    const ctx = canvas.getContext("2d");


    const ancho = canvas.clientWidth;

    const alto = canvas.clientHeight;

    if (ancho === 0 || alto === 0) {
        return;
    }


    const escala = window.devicePixelRatio || 1;

    canvas.width = Math.round(ancho * escala);

    canvas.height = Math.round(alto * escala);

    ctx.setTransform(escala, 0, 0, escala, 0, 0);


    ctx.clearRect(0, 0, ancho, alto);


    // Colores de los métodos.
    const series = [

        {
            nombre: "Euler",
            propiedad: propiedad === "velocidad" ? "eulerV" : "eulerY",
            color: "#ffb45c"
        },

        {
            nombre: "Heun",
            propiedad: propiedad === "velocidad" ? "heunV" : "heunY",
            color: "#4ba3ff"
        },

        {
            nombre: "RK4",
            propiedad: propiedad === "velocidad" ? "rk4V" : "rk4Y",
            color: "#4cdbad"
        },

        {
            nombre: "Analítica",
            propiedad: propiedad === "velocidad" ? "exactV" : "exactY",
            color: "#a58bff"
        }

    ];


    const margen = {
        izquierda: 60,
        derecha: 18,
        arriba: 20,
        abajo: 38
    };


    const anchoGrafica = ancho - margen.izquierda - margen.derecha;

    const altoGrafica = alto - margen.arriba - margen.abajo;


    const tiempoMaximo = filas[filas.length - 1].tiempo;


    const todosValores = [];

    series.forEach(serie => {

        filas.forEach(fila => {

            todosValores.push(fila[serie.propiedad]);

        });

    });


    let minimo = Math.min(0, ...todosValores);

    let maximo = Math.max(0, ...todosValores);


    if (minimo === maximo) {
        maximo = minimo + 1;
    }


    const margenVertical = (maximo - minimo) * 0.12;

    maximo += margenVertical;

    minimo -= margenVertical;


    // Posición horizontal.
    function posicionX(t) {

        return margen.izquierda
            + (t / tiempoMaximo) * anchoGrafica;

    }


    // Posición vertical.
    function posicionY(valor) {

        return margen.arriba
            + altoGrafica
            - ((valor - minimo) / (maximo - minimo)) * altoGrafica;

    }


    const estilo = getComputedStyle(document.body);

    const colorTexto = estilo.getPropertyValue("--text-secondary").trim();


    // Dibujar líneas horizontales y valores del eje Y.
    ctx.font = "10px Segoe UI";

    ctx.textAlign = "right";

    ctx.textBaseline = "middle";


    for (let i = 0; i <= 5; i++) {

        const valor = minimo + (maximo - minimo) * i / 5;

        const y = posicionY(valor);


        ctx.strokeStyle = "rgba(140,160,190,0.14)";

        ctx.lineWidth = 1;

        ctx.beginPath();

        ctx.moveTo(margen.izquierda, y);

        ctx.lineTo(ancho - margen.derecha, y);

        ctx.stroke();


        ctx.fillStyle = colorTexto;

        ctx.fillText(valor.toFixed(1), margen.izquierda - 10, y);

    }


    // Dibujar los valores del eje X.
    ctx.textAlign = "center";

    ctx.textBaseline = "top";


    for (let i = 0; i <= 5; i++) {

        const tiempo = tiempoMaximo * i / 5;

        const x = posicionX(tiempo);


        ctx.strokeStyle = "rgba(140,160,190,0.12)";

        ctx.beginPath();

        ctx.moveTo(x, margen.arriba);

        ctx.lineTo(x, alto - margen.abajo);

        ctx.stroke();


        ctx.fillStyle = colorTexto;

        ctx.fillText(tiempo.toFixed(1), x, alto - margen.abajo + 12);

    }


    // Dibujar los ejes principales.
    ctx.strokeStyle = colorTexto;

    ctx.lineWidth = 1;

    ctx.beginPath();

    ctx.moveTo(margen.izquierda, margen.arriba);

    ctx.lineTo(margen.izquierda, alto - margen.abajo);

    ctx.lineTo(ancho - margen.derecha, alto - margen.abajo);

    ctx.stroke();


    // Dibujar cada método.
    series.forEach(serie => {

        ctx.beginPath();

        ctx.strokeStyle = serie.color;

        ctx.lineWidth = 2.1;

        ctx.lineJoin = "round";

        ctx.lineCap = "round";


        filas.forEach((fila, indice) => {

            const x = posicionX(fila.tiempo);

            const y = posicionY(fila[serie.propiedad]);


            if (indice === 0) {

                ctx.moveTo(x, y);

            } else {

                ctx.lineTo(x, y);

            }

        });


        ctx.stroke();

    });


    // Etiquetas de los ejes.
    ctx.fillStyle = colorTexto;

    ctx.font = "10px Segoe UI";

    ctx.textAlign = "center";

    ctx.fillText("Tiempo (s)", ancho / 2, alto - 3);

}


// =========================================
// 16. MOSTRAR TODAS LAS GRÁFICAS
// =========================================

function mostrarGraficas(filas) {

    dibujarGrafica(
        "velocity-chart",
        "velocidad",
        filas
    );


    dibujarGrafica(
        "height-chart",
        "altura",
        filas
    );

}


// =========================================
// 17. EXPORTAR RESULTADOS A CSV
// =========================================

function exportarCSV() {

    if (resultadosActuales.length === 0) {

        mostrarMensaje("Primero debes ejecutar una simulación.");

        return;

    }


    const encabezado = [

        "Tiempo (s)",

        "Velocidad Euler (m/s)",
        "Altura Euler (m)",

        "Velocidad Heun (m/s)",
        "Altura Heun (m)",

        "Velocidad RK4 (m/s)",
        "Altura RK4 (m)",

        "Velocidad Analítica (m/s)",
        "Altura Analítica (m)"

    ];


    const filasCSV = resultadosActuales.map(fila => [

        fila.tiempo,
        fila.eulerV,
        fila.eulerY,
        fila.heunV,
        fila.heunY,
        fila.rk4V,
        fila.rk4Y,
        fila.exactV,
        fila.exactY

    ].map(valor => Number(valor).toFixed(8).replace(".", ",")));


    const contenido = "\uFEFF" + [
        encabezado.join(";"),
        ...filasCSV.map(fila => fila.join(";"))
    ].join("\n");


    const archivo = new Blob(
        [contenido],
        { type: "text/csv;charset=utf-8;" }
    );


    const enlace = document.createElement("a");

    const url = URL.createObjectURL(archivo);

    enlace.href = url;

    enlace.download = "resultados_euler_espacio.csv";


    document.body.appendChild(enlace);

    enlace.click();

    enlace.remove();

    URL.revokeObjectURL(url);

}


// =========================================
// 18. MENSAJES
// =========================================

function mostrarMensaje(texto, tipo = "error") {

    const mensaje = document.getElementById("message");

    mensaje.textContent = texto;

    mensaje.style.color = tipo === "error"
        ? "#ffb45c"
        : "#4cdbad";

}


// =========================================
// 19. EJECUTAR LA SIMULACIÓN COMPLETA
// =========================================

function calcularMision() {

    const parametros = leerParametros();


    // Comprobar los datos.
    const error = validarParametros(parametros);

    if (error !== "") {

        mostrarMensaje(error);

        return;

    }


    // Guardar la configuración.
    parametrosActuales = parametros;


    // Ejecutar los métodos.
    resultadosActuales = ejecutarSimulacion(parametros);


    // Mostrar los resultados.
    mostrarResultados(resultadosActuales);

    mostrarTabla(resultadosActuales);

    mostrarGraficas(resultadosActuales);


    mostrarMensaje(
        "Simulación completada correctamente.",
        "success"
    );

}


// =========================================
// 20. RESTABLECER LOS VALORES INICIALES
// =========================================

function restablecerValores() {

    document.getElementById("mass").value =
        valoresIniciales.masa;

    document.getElementById("thrust").value =
        valoresIniciales.empuje;

    document.getElementById("gravity").value =
        valoresIniciales.gravedad;

    document.getElementById("resistance").value =
        valoresIniciales.resistencia;

    document.getElementById("initial-velocity").value =
        valoresIniciales.velocidadInicial;

    document.getElementById("initial-height").value =
        valoresIniciales.alturaInicial;

    document.getElementById("total-time").value =
        valoresIniciales.tiempoTotal;

    document.getElementById("step-size").value =
        valoresIniciales.paso;


    calcularMision();

}


// =========================================
// 21. CAMBIAR ENTRE TEMA CLARO Y OSCURO
// =========================================

function configurarTema() {

    const boton = document.getElementById("theme-toggle");

    const temaGuardado = localStorage.getItem("tema-euler");


    if (temaGuardado === "claro") {

        document.body.classList.add("light-theme");

        boton.textContent = "☾";

    }


    boton.addEventListener("click", () => {

        document.body.classList.toggle("light-theme");


        const estaClaro = document.body.classList.contains("light-theme");


        localStorage.setItem(
            "tema-euler",
            estaClaro ? "claro" : "oscuro"
        );


        boton.textContent = estaClaro ? "☾" : "☼";


        // Volver a dibujar las gráficas para actualizar sus textos.
        if (resultadosActuales.length > 0) {

            mostrarGraficas(resultadosActuales);

        }

    });

}


// =========================================
// 22. CONECTAR LOS BOTONES
// =========================================

document.getElementById("parameters-form").addEventListener(
    "submit",
    function(evento) {

        evento.preventDefault();

        calcularMision();

    }
);


document.getElementById("reset-button").addEventListener(
    "click",
    restablecerValores
);


document.getElementById("export-button").addEventListener(
    "click",
    exportarCSV
);


window.addEventListener("resize", () => {

    if (resultadosActuales.length > 0) {

        mostrarGraficas(resultadosActuales);

    }

});


// =========================================
// 23. INICIAR LA APLICACIÓN
// =========================================

configurarTema();

calcularMision();
