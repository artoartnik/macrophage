// --- FUNCTiONS --- //
function polyArea(polygon){
  var A = 0;
  for(i=0; i<polygon.length-1; i++){
    A += polygon[i].x*polygon[i+1].y - polygon[i+1].x*polygon[i].y;
  }
  if( A < 0) A *= -1; 
  return A/2;
}

function unitVector(a, b){
  var vector = {};
  var D = dist(a.x, a.y, b.x, b.y);
  
  vector.x = (a.x-b.x)/D;
  vector.y = (a.y-b.y)/D;
  
  return vector;
}

function hookesLaw(a, b, D, k){
  var Acc = {};
  var X = dist(a.x, a.y, b.x, b.y);
  var F = (D-X)*k;
  var U = unitVector(a,b);
  
  Acc.ax = F*U.x;
  Acc.ay = F*U.y;
  
  return Acc;
}

//+ Jonas Raoni Soares Silva
//@ http://jsfromhell.com/math/is-point-in-poly [rev. #0]

function isPointInPoly(poly, pt){
    for(var c = false, i = -1, l = poly.length, j = l - 1; ++i < l; j = i)
        ((poly[i].y <= pt.y && pt.y < poly[j].y) || (poly[j].y <= pt.y && pt.y < poly[i].y))
        && (pt.x < (poly[j].x - poly[i].x) * (pt.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x)
        && (c = !c);
    return c;
}

// --- CLASSES --- //
function Vertex(x, y, px, py, vx, vy, ax, ay){
  this.x = x;
  this.y = y;
  this.px = px;
  this.py = py;
  this.vx = vx;
  this.vy = vy;
  this.ax = ax;
  this.ay = ay;
  this.links = [];
}

function Erytrocyte(x, y, vx, vy){
  this.x = x;
  this.y = y;
  this.vx = vx;
  this.vy = vy;
  
  this.move = function(){
    this.x += this.vx;
    this.y += this.vy;
  }
  
  this.draw = function(){
    noStroke();
    fill( 255, 150, 150 );
    ellipse(this.x,this.y,30,30);
  }
}

function Bacterium(x, y, vx, vy){
  this.x = x;
  this.y = y;
  this.vx = vx;
  this.vy = vy;
  this.ax = 0;
  this.ay = 0;
  this.dead = false;
  
  this.move = function(){
    this.x += this.vx;
    this.y += this.vy;
    
    this.vx += this.ax;
    this.vy += this.ay;
    
    this.ax = 0;
    this.ay = 0;
  }
  
  this.draw = function(){
    noStroke();
    if(!this.dead){
      fill( 16, 16, 16 );
      ellipse(this.x,this.y,10,10);
    } else {
      fill( 100, 100, 100 );
      ellipse(this.x,this.y,5,5);
    }   
  }
}

// Class for the macrophage
function Phage(nVertexes, x, y){
  this.membrane = [nVertexes];
  this.center = {"x": x, "y": y};
  this.moving = false;
  
  // generate regular n-polygon with linked vertices
  var degrees = TWO_PI/nVertexes;
  for(var i=0; i<nVertexes; i++){
    this.membrane[i] = new Vertex(x+100*sin(i*degrees),y+100*cos(i*degrees),x+100*sin(i*degrees),y+100*cos(i*degrees),0,0,0,0);
    
    // link with previous
    if(i > 0){
      this.membrane[i].links.push(this.membrane[i-1]);
      this.membrane[i-1].links.push(this.membrane[i]);
    }
    
    // link last with first
    if(i == (nVertexes-1)) {
      this.membrane[i].links.push(this.membrane[0]);
      this.membrane[0].links.push(this.membrane[i]);
    }
  }
    
  this.move = function(x, y, eatable){
      var H; // For use of hookes.
      
      // Find closest vertex in membrane
      var closestDist = Number.MAX_VALUE;
      var closestVertex = null;
      for(var i in this.membrane){
        var distToPoint = dist(this.membrane[i].x,this.membrane[i].y, x, y);
        if( distToPoint < closestDist ){
      
        closestVertex = this.membrane[i];
        closestDist = distToPoint;
        }
      }
      
      if(closestDist < 5){
        this.moving = false;
      }
      
      if(this.moving){
        closestVertex.ax = 10*(x - closestVertex.x)/closestDist;
        closestVertex.ay = 10*(y - closestVertex.y)/closestDist;
      }
      
      // Calculate forces
      for(var i in this.membrane){
        // Check for eatable objects and eat them
        for(var e in eatable){
          if(isPointInPoly(this.membrane, eatable[e])){
            eatable[e].dead = true;
          }                                        
        }
          
        for(var j in this.membrane[i].links){
          H = hookesLaw(this.membrane[i], this.membrane[i].links[j], 10, 0.05);
          this.membrane[i].ax += H.ax;
          this.membrane[i].ay += H.ay; 
        }
      
        H = hookesLaw(this.membrane[i], this.center, 80, 0.05);
        this.membrane[i].ax += H.ax;
        this.membrane[i].ay += H.ay;
      }
      
      // Move using Verlet integration
      for(var i in this.membrane){
        this.membrane[i].x = 2*this.membrane[i].x - this.membrane[i].px + this.membrane[i].ax*DT*DT;
        this.membrane[i].px = this.membrane[i].x;
        
        this.membrane[i].y = 2*this.membrane[i].y - this.membrane[i].py + this.membrane[i].ay*DT*DT;
        this.membrane[i].py = this.membrane[i].y;
        
         
        // Reset acceleration
        this.membrane[i].ax = 0;
        this.membrane[i].ay = 0;
     }
      
              
     // Move center
     var cX = 0
     var cY = 0;
     
     for(i in this.membrane){
      cX += this.membrane[i].x;
      cY += this.membrane[i].y;
     }
     
     this.center.x = cX/this.membrane.length;
     this.center.y = cY/this.membrane.length;
   }
   
   // Draw the macrophage                                   
   this.draw = function(){
      /*
      fill(255,255,255);
      beginShape();
      for(var i in this.membrane){
        vertex(this.membrane[i].x, this.membrane[i].y);
      }
      endShape();
      */
      
      fill( #6F50A1, 50 );
      noStroke();
      beginShape();
      
      // draw a curve between the points
      for(var i in this.membrane){
        curveVertex(this.membrane[i].x, this.membrane[i].y);
      }
      
      // close the curve
      for(i=0; i<3; i++){
      curveVertex(this.membrane[i].x, this.membrane[i].y);
      }
          
      endShape();
      
      fill( #6F50A1 );
      ellipse(this.center.x, this.center.y, 30, 30);
      
      /*
      for(var i in this.membrane){
        if(i == 0){
          fill(0,255,0);
        } else if( i == (this.membrane.length -1)){
          fill(0,0,255);
        } else {
          fill(255,0,0);
        }
        ellipse(this.membrane[i].x, this.membrane[i].y,5,5);
      }
      */
      
          
   }
}

// --- GLOBALS --- //
var DT = 1; //timestep

var PLAYER = new Phage(8, 200, 200);
var PMOVEX = 0
var PMOVEY = 0;

var CELLS = [];
var BACTERIA = [];
  
// --- PROCESSING.JS LOOP --- //

// Setup the Processing Canvas
void setup(){
  size( window.innerWidth, window.innerHeight );
  frameRate( 60 );
  colorMode(RGB, 255);
  background( #DDDDDD );
  
  for(var i=0; i<100; i++){
    CELLS.push(new Erytrocyte(random(0, window.innerWidth),random(0, window.innerHeight), random(-0.1,0.1), random(-0.1,0.1)));
  }
  
  for(var i=0; i<10; i++){
    BACTERIA.push(new Bacterium(random(0, window.innerWidth),random(0, window.innerHeight), random(-0.1,0.1), random(-0.1,0.1)));
  }
}

// Main draw loop
void draw(){           
  background( #DDDDDD );
  
    
  // Draw cells
  for(i in CELLS){
    CELLS[i].move();
    CELLS[i].draw();
  }
  
  // Draw bacteria
  for(i in BACTERIA){
    BACTERIA[i].move();
    BACTERIA[i].draw();
  }
  
  if( mousePressed ){
    PMOVEX = mouseX;
    PMOVEY = mouseY;
    PLAYER.moving = true;
  } else {
    PLAYER.moving = false;
  }
  
  PLAYER.move(PMOVEX, PMOVEY, BACTERIA);
  PLAYER.draw();
}