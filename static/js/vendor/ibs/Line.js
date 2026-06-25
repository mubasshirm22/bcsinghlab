//export {Line, Points, group_package, createLine};
//除了start end type point_num id之外的所有参数都会被传入style对象中
import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";
import { Points } from "./Points.js";
import IbsUtils from "./IbsUtils.js";
import { Rotater } from "./Rotater.js";

const defaultLine = {
  type: "straight", //straight, polyline, curve
  //如果position里点的数目超过两个，则以position为准，如果只有两个，则以pointsNum为准
  position: "0,0 100,0",
  pointsNum: 2,
  id: "",
  displayName: "",
  style: {
    isDash: false,
    left: "none",
    right: "none",
    color: "black",
    size: 1,
  },
  editable: true,
}
/**
 * create a line object
 * @param {{type:String, pointsPosition:String, id:String, displayName: String, style:{
 * isDash:Boolean[Number], left:String, right:String, color:String, size:Number}}} lineOption 
 * @param {*} parent 
 */
function Line(lineOption, parent) {
  this.parent = parent;
  let {
    option,
    config
  } = this.standarizeOptionAndConfig(lineOption);
  this.option = option;
  this.config = config;
  this.children = [];
  if ("site" !== parent.config.type) {
    this.nodesIndex = parent.children.length;
  }
  // this.nodesIndex = null;
  this.nodes = this.createNodes(); //该对象对应的svg节点
  // IbsUtils.showBorderByMousedown(this.nodes,this.config.containerNode);
  if ("site" !== parent.config.type) {
    this.draggable(this.option.editable);
  }
  
}
//根据对象创建path节点
/**
 * standarize the input line option
 * @param {{}} lineOption 
 */
Line.prototype.standarizeOptionAndConfig = function (lineOption, updateId = true) {
  let option = IbsUtils.supplementToStandardObj(defaultLine, lineOption);
  option.editable = this.parent.config.editable;
  let config = {};
  option.type = option.type.toLowerCase();
  if (updateId) {
   option.id = IbsUtils.getRandomId("Line"); 
  }
  option.displayName = option.displayName ? option.displayName : option.id;
  option.style.size = Math.max(1, parseInt(option.style.size));
  // let CoorArray = [];
  /* for (let i in this.option.coordiante) {
    CoorArray.push([this.option.coordinate[i].x,this.option.coordinate[i],y].join(","))
  } */
  if (typeof option.position == "object") {
    let coor = [];
    for (let key in option.position) {
      let stringCoor = [option.position[key].x, option.position[key].y].join(",");
      coor.push(stringCoor);
    }
    option.position = coor.join(" ")
  }
  let CoorArray = option.position.split(" "); //["1,2","3,4"]
  let pointArray = []; //[[1,2],[3,4]]
  for (let i =0; i < CoorArray.length; i++) {
    let XandY = CoorArray[i].split(",");
    pointArray.push([parseFloat(XandY[0]), parseFloat(XandY[1])]);
  }
  let pathArray = [];
  if ("straight" === option.type) {
    pathArray[0] = pointArray[0];
    pathArray[1] = pointArray.pop();
  } else {
    if (CoorArray.length == option.pointsNum) {
      pathArray = pointArray;
    } else {
      let len = CoorArray.length;
      for (let i = 0; i < option.pointsNum; i++) {
        let x = pointArray[0][0] + (pointArray[len - 1][0] - pointArray[0][0]) / (option.pointsNum - 1) * i
        let y = pointArray[0][1] + (pointArray[len - 1][1] - pointArray[0][1]) / (option.pointsNum - 1) * i
        pathArray.push([x, y]);
      }
    }
  }
  
  let newPosition = "";
  pathArray.forEach(coor => {
    newPosition += `${coor[0]},${coor[1]} `;
  });
  option.position = newPosition.trim();
  
  pathArray[0].unshift("M");
  if (option.type == "polyline" || option.type == "straight") {
    pathArray[1].unshift("L");
  } else if (option.type == "curve" && (parseInt(option.pointsNum) - 3) % 2 == 0) {
    pathArray[1].unshift("Q");
  } else if (option.type == "curve" && (parseInt(option.pointsNum) - 4) % 3 == 0) {
    pathArray[1].unshift("C");
  } else {
    throw ("pointsNum or line type error.")
  }
  option.pointsNum = pathArray.length;
  config.pathArray = pathArray;
  config.type = "line";
  config.rotate = 0;
  return {
    option,
    config
  };
}
/**
 * create nodes of the line object
 */
Line.prototype.createNodes = function () {
  let styleOption = IbsUtils.standarizeStyleOption(this.option.style);
  let lineAttr = {
    d: this.produceD(),
    id: this.option.id
  }
  let line = IbsUtils.createSvgElement("path", lineAttr);
  for (let key in styleOption) {
    line.style[key] = styleOption[key];
  }
  line.style.fill = "none";
  line.setAttribute("id", this.option.id);
  this.config.lineNode = line;
  let group = IbsUtils.createSvgElement("g", {
    id: `${this.option.id}Group`
  });
  group.appendChild(line);
  if (this.option.editable) {
    let outline = this.createOutLine(line);
    let container = this.createControlPoint();
    this.config.containerNode = container;
    container.setAttribute("visibility", "hidden");
    group.appendChild(outline);
    group.appendChild(container);
    if ("root" === this.parent.config.type) {
      IbsUtils.showBorderByMousedown(group, this.config.containerNode);
      IbsUtils.showBorderAlone(this.config.containerNode);
      this.parent.config.selected = this;
      group.addEventListener("mousedown", () => {
        this.parent.config.selected = this;
        this.parent.config.isHideAll = false;
      })
    }
  }
  if (this.option.style.right != "none" || this.option.style.left != "none") {
    let defs = this.createMarker(line);
    this.config.defs = defs;
    group.appendChild(defs);
  }


  return group;
};
/**
 * create points used to change the line
 */
Line.prototype.createControlPoint = function () {
  for (let i = 0; i < this.config.pathArray.length; i++) {
    let XandY = this.config.pathArray[i].slice(-2);
    let pointOption = {
      cx: XandY[0],
      cy: XandY[1],
      r: 4,
      color: "#ffa500",
      index: i,
      editPoint: this.option.editable,
    }
    let point = new Points(pointOption, this);
    this.children.push(point);
  }
  let container = IbsUtils.createSvgElement("g", {
    class: "nodeBorder"
  });
  let start = this.config.pathArray[0].slice(-2);
  let end = this.config.pathArray.slice(-1)[0].slice(-2);
  let rCenterX = (start[0] + end[0]) / 2;
  let rCenterY = (start[1] + end[1]) / 2;
  let rotater = new Rotater({
    x: rCenterX,
    y: rCenterY,
    rotate: 0
  }, this);
  this.children.push(rotater);
  for (let i = 0; i < this.children.length; i++) {
    container.appendChild(this.children[i].nodes);
  }
  container.appendChild(rotater.config.cloneNode);
  return container
}
/**
 * create a outline which making easy to select the line 
 * @param {HTMLElement} lineNode 
 */
Line.prototype.createOutLine = function (lineNode) {
  let outline = lineNode.cloneNode(); //以克隆的形式生成外轮廓节点
  outline.setAttribute("id", `${this.option.id}Outline`);
  outline.style.cssText = "stroke:black;stroke-width:15;stroke-opacity:0;fill:none"
  this.config.outlineNode = outline;
  return outline;
}
/**
 * make the line draggable
 * @param {Boolean} editable 
 */
Line.prototype.draggable = function (editable = true) {
  let group_drag = SVG(this.nodes);
  if (editable) {
    group_drag.draggable();
    let dArray;
    // console.log(dArray);
    group_drag.on(`beforedrag.${this.option.id}Group`, (e) => {     
      dArray = SVG(this.config.lineNode).array();
      this.parent.config.isUndo = false;
      let operation = {
        target: this,
        cmd: "drag",
        args: JSON.stringify(this.option),
        nodesIndex: this.nodesIndex
      }
      this.parent.config.undoStack.push(operation);
    });
    group_drag.on(`dragmove.${this.option.id}Group`, (e) => {
      let {
        handler,
        box
      } = e.detail;
      e.preventDefault();
      handler.move(box.x, box.y);
      this.updatePathArray(0, dArray[0][1], dArray[0][2])
      if (dArray.length == 2) { // straight/curve
        for (let i = 1; i < this.config.pathArray.length; i++) {
          this.updatePathArray(i, dArray[1][2 * i - 1], dArray[1][2 * i]);
        }
      } else { // polyline
        for (let i = 1; i < dArray.length; i++) {
          this.updatePathArray(i, dArray[i][1], dArray[i][2]);
        }
      }
    })
    group_drag.on(`dragend.${this.option.id}Group`, () => {
      for (let i = 0; i < this.children.length; i++) {
        //只有rotater有config属性，points是没有的
        if (this.children[i].config) {
          let start = this.config.pathArray[0].slice(-2);
          let end = this.config.pathArray.slice(-1)[0].slice(-2);
          let rCenterX = (start[0] + end[0]) / 2;
          let rCenterY = (start[1] + end[1]) / 2;
          this.children[i].updateCoordinate(rCenterX, rCenterY);
          continue;
        }
        this.children[i].updateCoordinate(true);
      }
    }) 
  } else {
    group_drag.draggable(false);
  }
}
/**
 * update the node of line node
 * @param {{type:String, pointsPosition:String, id:String, displayName: String, style:{ isDash:Boolean[Number], left:String, right:String, color:String, size:Number}} lineOption
 */
Line.prototype.update = function (lineOption) {
  let previousOption = JSON.stringify(this.option);
  let {
    option,
    config
  } = this.standarizeOptionAndConfig(lineOption, false);
  // let pathArray = [];
  this.option = option;
  this.config = config;
  this.children = [];
  let updateNode = this.createNodes();
  this.parent.nodes.replaceChild(updateNode, this.nodes);
  this.nodes = updateNode;
  IbsUtils.showBorderByMousedown(this.nodes, this.config.containerNode)
  this.draggable();
  this.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "update",
    args: previousOption,
    nodesIndex: this.nodesIndex
  }
  this.parent.config.undoStack.push(operation);
  return this;
}
Line.prototype.delete = function () {
  IbsUtils.deleteObj(this);
  this.parent.config.isUndo = false;
  // this.parent.config.elementsCount.line -= 1;
  let operation = {
    target: this,
    cmd: "delete",
    args: JSON.stringify(this.option),
    nodesIndex: this.nodesIndex
  }
  this.parent.config.undoStack.push(operation);
  return true;
}
/**
 * enable change the control points coordinate based on given index
 * @param {{}} obj the object to change
 * @param {Number} x new x coordinate
 * @param {Number} y new y coordinate
 */
Line.prototype.updatePathArray = function (i, x, y, reshape = false) {
  x = Math.round(100 * x) / 100;
  y = Math.round(100 * y) / 100;
  this.config.pathArray[i].splice(-2, 1, x);
  this.config.pathArray[i].splice(-1, 1, y);
  //获取path对象的d属性矩阵
  if (reshape) {
    // console.log(SVG(this.config.lineNode).PathArray());
    var dArray = SVG(this.config.lineNode).array();
    if (this.config.outlineNode) {
      var outlineArray = SVG(this.config.outlineNode).array();
    }
    //在点拖动时更新d属性矩阵中对应的坐标
    switch (parseFloat(i)) {
      case 0:
        dArray[0][1] = x;
        dArray[0][2] = y;
        if (this.config.outlineNode) {
          outlineArray[0][1] = x;
          outlineArray[0][2] = y;
        }
        break;
      default:
        let flag = this.option.type == "polyline" || this.option.type == "straight";
        //console.log(flag);
        if (flag) {
          dArray[i][1] = x;
          dArray[i][2] = y;
          if (this.config.outlineNode) {
            outlineArray[i][1] = x;
            outlineArray[i][2] = y;
          }

        } else if (this.option.type == "curve") {
          dArray[1][2 * i - 1] = x;
          dArray[1][2 * i] = y;
          if (this.config.outlineNode) {
            outlineArray[1][2 * i - 1] = x;
            outlineArray[1][2 * i] = y;
          }
        }
    }
  }
  let positionArray = [];
  for (let i = 0; i < this.config.pathArray.length; i++) {
    let coorArray = this.config.pathArray[i].slice(-2);
    positionArray[i] = coorArray.join(",");
  }
  this.option.position = positionArray.join(" ");
}
/**
 * reshape the line while dragging the control points
 * @param {Number} index 
 * @param {Number} x 
 * @param {Number} y 
 * @param {Boolean} reshape 
 */
Line.prototype.reshape = function (index, x, y, reshape) {
  this.updatePathArray(index, x, y, reshape);
  let d = this.produceD();
  this.config.lineNode.setAttribute("d", d);
  if (this.option.editable) {
    this.config.outlineNode.setAttribute("d", d);
    if (index == 0 || index == this.config.pathArray.length - 1) {
      let start = this.config.pathArray[0].slice(-2);
      let end = this.config.pathArray.slice(-1)[0].slice(-2);
      let rCenterX = (start[0] + end[0]) / 2;
      let rCenterY = (start[1] + end[1]) / 2;
      this.children[this.children.length - 1].updateCoordinate(rCenterX, rCenterY);
    }
  }
}
/**
 * produce the attribute d based on the this.control_points
 */
Line.prototype.produceD = function () {
  var d = "";
  for (let i = 0; i < this.config.pathArray.length; i++) {
    d += this.config.pathArray[i].join(" ") + " ";
  }
  return d;
}
/**
 * create the markers on the one side of line based on the this.left or this.right
 * @param {String} position the marker position, left(start) or right(end) 
 */
Line.prototype.createMarker = function (lineNode) {
  let defs = document.createElementNS("http://www.w3.org/2000/svg", "defs");
  defs.setAttribute("id", `${this.option.id}Marker`);
  let markerOption = {
    id: `${this.option.id}Marker`,
    markerWidth: 100,
    markerHeight: 100,
    viewBox: "-10 -10 100 100",
    orient: "auto",
    refX: 0,
    refY: 0,
    markerUnits: "userSpaceOnUse"
  }
  if (this.option.style.left != "none") {
    let leftMarker = this.createOneSideMarker("left", markerOption);
    defs.appendChild(leftMarker);
    lineNode.style.markerStart = `url("#leftMarker${this.option.id}")`;
  }
  if (this.option.style.right != "none") {
    let rightMarker = this.createOneSideMarker("right", markerOption);
    defs.appendChild(rightMarker);
    lineNode.style.markerEnd = `url("#rightMarker${this.option.id}")`;
  }
  return defs;
}
/**
 *  create an arrow marker on the given position
 * @param {HTMLElement} markerNode the arrow will be the child of this node
 * @param {String} position define the arrow position on the line
 * @param {String} points define the arrow shape
 */
Line.prototype.createArrow = function (markerNode, position, points = "0 -2 3.26 0 0 2") { //创建箭头Marker
  var arrow = document.createElementNS("http://www.w3.org/2000/svg", "polygon"); //用三角形作为箭头
  //arrow属性设置
  arrow.setAttribute("stroke-width", this.option.style.size);
  arrow.setAttribute("points", points); //箭头形状设计，
  if (position == "left") { //如果是路径起点的标记，则将图形旋转180度
    arrow.setAttribute("transform", "rotate(180)");
  }
  markerNode.appendChild(arrow);
}
/**
 * create an scale marker on the given position
 * @param {HTMLElement} markerNode the arrow will be the child of this node
 * @param {String} position define the scale position on the line
 * @param {Number} x the horizon position of scale marker
 */
Line.prototype.createScale = function (markerNode, position, x = 0) {
  var scale = document.createElementNS("http://www.w3.org/2000/svg", "rect"); //以矩形作为标尺图形
  var h = this.option.style.size <= 5 ? 6 : this.option.style.size * 1.5;
  scale.setAttribute("x", x); //标尺左上角的x坐标
  scale.setAttribute("y", h * (-0.5)); //标尺左上角的y坐标，向上偏移高度的一半能够保证在标尺和线的垂直中心重合
  scale.setAttribute("width", this.option.style.size);
  scale.setAttribute("height", h);
  scale.style.cssText = `fill:${this.option.style.color}`;
  if (position == "left") { //如果是路径起点的标记，则旋转180度
    scale.setAttribute("transform", "rotate(180)");
  }
  markerNode.appendChild(scale);
}
Line.prototype.createOneSideMarker = function (side, markerOption) {
  let marker = IbsUtils.createSvgElement("marker", markerOption);
  marker.setAttribute("id", `${side}Marker${this.option.id}`);
  marker.style.cssText = `stroke-width:0.1;stroke:${this.option.style.color};fill:${this.option.style.color}`;
  switch (this.option.style[side].toLowerCase()) {
    case "arrow":
      this.createArrow(marker, side);
      break;
    case "scale":
      this.createScale(marker, side);
      break;
    case "arrowscale":
      this.createArrow(marker, side);
      this.createScale(marker, side, 3.26);
      break;
  }
  return marker;
}

Line.prototype.rotate = function (angel) {
  let rotater = this.children[this.children.length - 1];
  let rotate = angel - this.config.rotate;
  this.config.rotate = angel;
  let cx = rotater.option.x;
  let cy = rotater.option.y;
  for (let i = 0; i < this.config.pathArray.length; i++) {
    let coor = this.config.pathArray[i].slice(-2);
    let result = IbsUtils.calculatePointAfterRotate(cx, cy, coor[0], coor[1], rotate);
    this.reshape(i,result[0],result[1],true);
    this.children[i].updateCoordinate(false,coor[0],coor[1])
  }
}

export default Line
