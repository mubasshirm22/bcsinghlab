const pointNumTable = {
  pentagon: 5,
  hexagon: 6,
  octagon: 8
}
const initAngelTable = {
  pentagon: 0,
  hexagon: parseFloat(Math.PI) / 6,
  octagon: 0
}

/**
 * create triangle points coordiante
 */
function Triangle(control_points) {
  let points = [];
  points[0] = control_points[0];
  points[1] = control_points[3];
  points[2] = control_points[6];
  return {
    points
  };
}
/**
 * create points coordiante of pentagon, hexagon, octagon
 */
function Polygon(control_points) {
  let points = [];
  let pointNum = pointNumTable[this.option.shape];
  var angleStep = (2 * parseFloat(Math.PI)) / pointNum;
  let initAngel = initAngelTable[this.option.shape];
  var x;
  var y;
  var x0 = this.config.nodePosition.cx;
  var y0 = this.config.nodePosition.cy;
  var x_flag = x0 - control_points[0][0] > 0 ? 1 : -1;
  for (let i = 0; i < pointNum; i++) {
    var angel = initAngel + angleStep * i;
    x = x0 + this.option.width * 0.5 * Math.sin(angel) * x_flag;
    y = y0 - this.option.height * 0.5 * Math.cos(angel);
    points[i] = [x, y];
  }
  return {
    points
  };
}

/**
 * create points coordiante of rhombus
 */
function Rhombus(control_points) {
  let points = [];
  let order = [1, 3, 5, 7];
  for (let i = 0; i < order.length; i++) {
    points[i] = control_points[order[i]]
  }
  return {
    points
  };
}
/**
 * create points coordiante of parallelogram
 */
function Parallelogram(control_points) {
  let points = new Array(4);
  for (let i = 0;i<4;i++) {
    points[i] = [0,0];
  }
  // console.log(points);
  var x0 = this.config.nodePosition.cx;
  var x_flag = (x0 - control_points[0][0]) > 0 ? 1 : (-1);
  points[0] = control_points[0];
  points[1][0] = control_points[0][0] + 0.75 * this.option.width ;
  points[1][1] = control_points[0][1];
  points[2] = control_points[4];
  points[3][0] = control_points[6][0] + 0.25 * this.option.width ;
  points[3][1] = control_points[6][1];
  return {
    points
  };
}
/**
 * create points coordiante of trapezium
 */
function Trapezium(control_points) {
  let points = new Array(4);
  for (let i = 0;i<4;i++) {
    points[i] = [0,0];
  }
  var x0 = this.config.nodePosition.cx;
  var x_flag = (x0 - control_points[0][0]) > 0 ? 1 : (-1);
  points[0] = control_points[0];
  points[1] = control_points[2];
  points[2][0] = control_points[4][0] - this.option.width * 0.25 * x_flag;
  points[2][1] = control_points[6][1];
  points[3][0] = control_points[6][0] + this.option.width * 0.25 * x_flag;
  points[3][1] = control_points[6][1];
  return {
    points
  }
}
/**
 * create points coordiante of arrow
 */
function Arrow(control_points) {
  let points = new Array(7);
  for (let i = 0;i<7;i++) {
    points[i] = [0,0];
  }
  var x0 = this.config.nodePosition.cx;
  var x_flag = (x0 - control_points[0][0]) > 0 ? 1 : (-1);
  points[0][0] = control_points[7][0];
  points[0][1] = control_points[7][1] - this.option.height * 0.5 / 1.618 * 0.618;
  points[1][0] = control_points[7][0] + this.option.width / 1.618 * x_flag
  points[1][1] = control_points[7][1] - this.option.height * 0.5 / 1.618 * 0.618
  points[2][0] = control_points[7][0] + this.option.width / 1.618 * x_flag;
  points[2][1] = control_points[0][1]
  points[3] = control_points[3]
  points[4][0] = control_points[7][0] + this.option.width / 1.618 * x_flag
  points[4][1] = control_points[6][1]
  points[5][0] = control_points[7][0] + this.option.width / 1.618 * x_flag
  points[5][1] = control_points[7][1] + this.option.height * 0.5 / 1.618 * 0.618
  points[6][0] = control_points[7][0]
  points[6][1] = control_points[7][1] + this.option.height * 0.5 / 1.618 * 0.618
  return {points};
}
/**
 * create points coordiante of dovetail arrow
 */
function DovetailArrow(control_points) {
  let points = [];
  let order = [0, 1, 3, 5, 6];
  for (let i = 0; i < order.length; i++) {
    points[i] = control_points[order[i]]
  }
  points.push([this.config.nodePosition.cx,this.config.nodePosition.cy]);
  return {points};
}
/**
 * use to produce rectangle option 
 */
function Rect() {
  let contentOption = {};
  if (this.config.type == "marker") {
    contentOption = {
      x: this.config.nodePosition.cx - this.option.width * 0.5,
      y: this.config.nodePosition.cy - this.option.height * 0.5,
      width: this.option.width,
      height: this.option.height,
    }
  } else {
    let side = this.config.sidePoints;
    contentOption = {
      points:[["M",side[0][0],side[0][1]],["L",side[2][0],side[2][1]],["L",side[4][0],side[4][1]],["L",side[6][0],side[6][1]],["L",side[0][0],side[0][1]]]
    }
  }
  return contentOption;
}
function RoundRect() {
  if (this.config.type == "marker") {
    var contentOption = {
      x: this.config.nodePosition.cx - this.option.width * 0.5,
      y: this.config.nodePosition.cy - this.option.height * 0.5,
      width: this.option.width,
      height: this.option.height,
      rx: 3
    }
  } else {
    let side = this.config.sidePoints;
    let r= 3;
    let ds = 3 * Math.sin(this.option.style.rotate/180 * Math.PI);
    let dc = 3 * Math.cos(this.option.style.rotate/ 180 *Math.PI);
    let a = 0;
    // let index = [[0,1,1],[0,1,1],[0,0,0]]
    let points = [["M", side[0][0]+dc,side[0][1]+ds],["L",side[2][0]-dc,side[2][1]-ds],["A",r,r,a,0,1,side[2][0]-ds,side[2][1]+dc],["L",side[4][0]+ds,side[4][1]-dc],["A",r,r,a,0,1,side[4][0] - dc,side[4][1] - ds],["L",side[6][0]+dc,side[6][1] + ds],["A",r,r,a,0,1,side[6][0] + ds,side[6][1]-dc],["L",side[0][0] - ds,side[0][1]+dc],["A",r,r,a,0,1,side[0][0]+dc,side[0][1]+ds]]
    // let points = [];
    // for (let i = 0;i<7;i = i +2) {
    //   points[i/2] = this.config.sidePoints[i]
    // }
    var contentOption = {
      points,
    }
  }
  return contentOption;
}
/**
 * use to produce circle option
 */
function Circle() {
  if (this.config.type == "marker" ) {
    var contentOption = {
      cx: this.config.nodePosition.cx,
      cy: this.config.nodePosition.cy,
      rx: this.option.width / 2,
      ry: this.option.height / 2
    } 
  } else if (this.config.type == "domain" || this.config.type== "site") {
    let side = this.config.sidePoints;
    let points = [["M",side[7][0],side[7][1]],["A",this.option.width/2,this.option.height/2,this.option.style.rotate,1,1,side[3][0],side[3][1]],["A",this.option.width/2,this.option.height/2,this.option.style.rotate,1,1,side[7][0],side[7][1]]];
    var contentOption = {
      points
    }
  }
  
  return contentOption;
}

function JustAPoint() {
  var contentOption = {
    cx: this.config.nodePosition.cx,
    cy: this.config.nodePosition.cy,
    r:3,
    // visibility: "hidden",
    opacity: 0
  }
  return contentOption;
}

function Cylinder() {
  let points = this.config.sidePoints;
  var contentOption = {
    // d:`M${points[0][0]} ${points[0][1]} A 80 ${this.option.height} 0 0 1 ${points[6][0]} ${points[6][1]} L ${points[4][0]} ${points[4][1]} A 80 ${this.option.height} 0 0 0 ${points[2][0]} ${points[2][1]} L ${points[0][0]} ${points[0][1]} A 80 ${this.option.height} 0 0 0 ${points[6][0]} ${points[6][1]}`
    // d:`M${points[7][0] - 6} ${points[7][1] +2.5} A 6 ${this.option.height/2} 0 1 0 ${points[7][0] - 6} ${points[7][1] -2.5} h 4 a 5 4 0 0 1 0 5 h -4 M ${points[0][0]} ${points[0][1]} L ${points[2][0]} ${points[2][1]} A 6 ${this.option.height/2} 0 0 1 ${points[4][0]} ${points[4][1]} L ${points[6][0]} ${points[6][1]} A 6 ${this.option.height/2} 0 0 0 ${points[0][0]} ${points[0][1]}`
    points:[["M",points[7][0],points[7][1] +2.5], ["A",6, this.option.height/2,0,1,0,points[7][0],points[7][1] -2.5],["H",points[7][0]+4],["A",5,4,0,0,1,points[7][0]+4,points[7][1] + 2.5],["H",points[7][0]],["M",points[0][0] + 6,points[0][1]],["L",points[2][0] -6,points[2][1]],["A",6,this.option.height/2,0,0,1,points[4][0] -6,points[4][1]],["L",points[6][0] + 6,points[6][1]],["A",6,this.option.height/2,0,0,0,points[0][0] + 6,points[0][1]]],
    // d:`M${points[7][0]} ${points[7][1] +2.5} A 6 ${this.option.height/2} 0 1 0 ${points[7][0]} ${points[7][1] -2.5} h 4 a 5 4 0 0 1 0 5 h -4 M ${points[0][0] + 6} ${points[0][1]} L ${points[2][0] -6} ${points[2][1]} A 6 ${this.option.height/2} 0 0 1 ${points[4][0] -6} ${points[4][1]} L ${points[6][0] + 6} ${points[6][1]} A 6 ${this.option.height/2} 0 0 0 ${points[0][0] + 6} ${points[0][1]}`  
  }
  return contentOption;
}

const methodsTable = {
  triangle: Triangle,
  pentagon: Polygon,
  hexagon: Polygon,
  octagon: Polygon,
  rhombus: Rhombus,
  parallelogram: Parallelogram,
  trapezium: Trapezium, //
  arrow: Arrow, //
  dovetailArrow: DovetailArrow,
  rect: Rect,
  roundRect: RoundRect,
  circle: Circle,
  "none":JustAPoint,
  cylinder:Cylinder
}

export {
  Triangle,
  Polygon,
  Rhombus,
  Parallelogram,
  Trapezium,
  Arrow,
  DovetailArrow,
  Rect,
  Circle,
  RoundRect,
  methodsTable,
  Cylinder
}
