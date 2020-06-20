var params = {
    "L": 2, "H":1, "mu":0, "k":0
} ;

function updateValText(slider_id, value){
    this_name = slider_id.replace("_slider", "");
    text = document.getElementById(this_name);
    text.innerHTML = value;
}

function parseSlider(name){
    // Take a name, parse slider with given id to update text and parameter on update
    slider = document.getElementById(name+"_slider");
    updateValText(slider.id, slider.value);
    slider.oninput = function() {
        updateValText(this.id, this.value);
        params[this_name] = parseFloat(this.value);
        plotBrach()
    }
    //slider.oninput();
}

// Initial setup
var canvas = document.getElementById("simAll");
const ctx = canvas.getContext("2d");
var width = canvas.clientWidth;
var height = canvas.clientHeight;

function setUpCanvas() {
    // set up canvas

    width = canvas.clientWidth;
    height = canvas.clientHeight;
    canvas.width = width;
    canvas.height = height;
    ctx.fillStyle = 'rgb(255, 255, 255)';
    ctx.fillRect(0, 0, width, height);

    canvas_x0 = 0.01; // as a frac of width
    canvas_y0 = 0.01;
    canvas_xwidth = 0.98;
    canvas_yheight = 0.98;

    var max_H = 3;
    var max_L = 3;

    xmin = 0;
    ymin = 0;
    xmax = max_L+0.01;
    ymax = max_H+0.01;

}

function drawCurve(X, Y, c="#4CC984", lw=5) {
// Given an array of X and Y values, plot a curve of the points + circles on either end
// color default is blue
    var x, y;
    ctx.beginPath();
    for (i=0; i<X.length; i++) {
        x = X[i];
        y = Y[i];
        ctx.lineTo(x, y);
    }
    ctx.lineWidth = lw;
    ctx.strokeStyle = c;
    ctx.setLineDash([]);
    ctx.stroke();

    // Draw circle for end points
    for (i of [0, X.length-1]) {
        ctx.beginPath();
        ctx.arc(X[i], Y[i], 5, 0, 2 * Math.PI);
        ctx.fillStyle = c;
        ctx.fill()
    }

}

function drawError(L, H, lw=5){
    // Draw error to show no valid curve
    ctx.beginPath();
    var A = coordTransform([0, 0])
    var B = coordTransform([L,H])
    ctx.lineTo(A[0], A[1])
    ctx.lineTo(B[0], B[1]);
    ctx.lineWidth = lw;
    ctx.strokeStyle = "red";
    // ctx.setLineDash([width/20,width/20]);
    ctx.stroke();

}

function coordTransform(coord){
    // convert x y from equation coordinates to canvas coordinates
    x = coord[0];
    y = coord[1];
    var x_glob = (canvas_x0 + ((x-xmin)/(xmax-xmin)) * canvas_xwidth) * width;
    var y_glob = (canvas_y0 + ((y-ymin)/(ymax-ymin)) * canvas_yheight) * height;
    return [x_glob, y_glob]
}

function reverseCoordTransform(coord){
    // convert x y from canvas coordinates to equation coordinates
    x = coord[0];
    y = coord[1];
    var x_loc = (((x / width) - canvas_x0) * (xmax-xmin) / canvas_xwidth) + xmin;
    var y_loc = (((y / height) - canvas_y0) * (ymax-ymin) / canvas_yheight) + ymin;
    return [x_loc, y_loc]
}

function brach(theta, mu=0, A=0.5){
    var x = A * (theta - Math.sin(theta) + mu * (1-Math.cos(theta)));
    var y = A * (1- Math.cos(theta) + mu * (theta+Math.sin(theta)));

    return [x, y]
}

function getPhi(L, H, mu=0, step=0.25, tol=1e-4, maxit=500,
) {
    // For a given H, find the value of theta which means brach() goes through the point (L,H)
    // Return this as phi

    // Function to find root of
    function f(theta) {
        coord = brach(theta, mu); // independent of A
        return coord[1]/coord[0] - H/L;
    }

    // Solve in 2 stages: increment by step at a time until root is crossed, then
    // Interval bisection. This prevents skipping of the first root
    // Solve by interval bisection
    var phi_low = 1e-3;
    var phi_high = step;
    for(i=0; phi_high<10; i++) {
        phi_high += step;
        if (Math.sign(f(phi_low)) !== Math.sign(f(phi_high))){
            phi_low = phi_high - step;
            break
        }
    }
    // If no range found, use full 0-10 range for interval bisection
    if (phi_high === 10){
        phi_low = 0
    }

    var mid;
    for (i=0; i<maxit; i++){
        mid = ( phi_low + phi_high ) / 2;
        if (Math.sign(f(phi_high)) === Math.sign(f(mid)))
        {
            phi_high = mid;
        }
        else {
            phi_low = mid;
        }

        // End if reached sufficient accuracy
        if (Math.abs(f(mid)) <= tol) {
            return mid;
        }
    }
    return mid;
}

function getA(phi, H, mu=0){

    var A = H / (1- Math.cos(phi) + mu * (phi+Math.sin(phi)));
    return A
}

function plotBrach() {
    ctx.clearRect(0, 0, width, height);
    var X = [];
    var Y = [];
    var mu = params["mu"];
    var L = params["L"];
    var H = params["H"];

    var phi = getPhi(L, H, mu);
    if (phi<1e-2){
        phi = getPhi(L,H,mu, 0.1); // attempt again with finer step
    }

    var A = getA(phi, H, mu);

    document.getElementById("phi").innerHTML = phi.toString();
    document.getElementById("A").innerHTML = A.toString();

    if (phi>1e-2) {

        for (theta = 0; theta < phi; theta += 0.01) {
            coords = coordTransform(brach(theta, mu, A));
            X.push(coords[0]);
            Y.push(coords[1]);
        }
        drawCurve(X, Y);
    }
    else{
        drawError(L, H)
    }
}


var mouse_isdown = false;
// Add event listener for `click` events.
canvas.addEventListener('mousedown', function(event) {
    mouse_isdown = true;
}, false);

canvas.addEventListener('mousemove', function(event) {
    if (mouse_isdown){
    var x = event.pageX - canvas.offsetLeft,
    y = event.pageY - canvas.offsetTop;

    coords = reverseCoordTransform([x, y]);

    params["L"] = coords[0];
    params["H"] = coords[1];
    plotBrach()

    }

}, false);

canvas.addEventListener('mouseup', function(event) {
    mouse_isdown = false;

}, false);

// On window resize, resize canvas
window.addEventListener('resize', setUpCanvas);

setUpCanvas();

slider_names = ["mu", "k"];
for (name of slider_names) {
    parseSlider(name)
}

plotBrach();

