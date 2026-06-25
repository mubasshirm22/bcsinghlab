import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";
import { Points } from "./Points.js";
import IbsUtils from "./IbsUtils.js";
import { Rotater } from "./Rotater.js";

/**
 * @param {{width:Number, height:Number,id:String, position:{cx:Number, cy:Number}}} containerOption
 * @param {IbsCharts} parent
 */
function Container(containerOption, parent) {
  this.parent = parent;
  this.children = []; //point * 8, rotater
  let { option, config } = this.standarizeOptionAndConfig(containerOption);
  this.count = 0;
  this.option = option;
  this.config = config;
  this.nodes = this.createNode();
}
Container.prototype.standarizeOptionAndConfig = function (containerOption) {
  let option = containerOption;
  let config = {};
  let { sidePoints, centerPoint } = IbsUtils.calRectEightPoints(
    option.position.cx,
    option.position.cy,
    option.width,
    option.height,
    this.parent.option.style.rotate
  );
  config.sidePoints = sidePoints;
  config.centerPoint = centerPoint;
  config.nodesTransform = [0, 0, 0];
  config.type = "container";
  config.coorCenter = {
    x: option.position.cx,
    y: option.position.cy,
  };
  return {
    option,
    config,
  };
};
/**
 * create nodes and set its attribute and style
 */
Container.prototype.createNode = function () {
  let rectOption = {
    points: this.producePointsAttribute(),
    // id:this.parent.option.id + "Container",
    style: "fill:none;fill-opacity:0;stroke:black",
  };
  let rect = IbsUtils.createSvgElement("polygon", rectOption);
  this.config.rectNode = rect;
  rect.setAttribute("visibility", "hidden");
  let group = IbsUtils.createSvgElement("g", {
    id: this.parent.option.id + "ContainerGroup",
    class: "nodeBorder",
  });
  group.appendChild(rect);
  this.createControlPoints();
  if (this.parent.option.shape != "cylinder") {
    var rotater = this.createRotater();
    this.children.push(rotater);
  }
  for (let i = 0; i < this.children.length; i++) {
    group.appendChild(this.children[i].nodes);
    if (parseFloat(i) % 2 != 0) {
      this.children[i].nodes.setAttribute("visibility", "hidden");
    }
  }
  return group;
};
/**
 * create 8 control points for the container.
 */
Container.prototype.createControlPoints = function () {
  for (let i = 0; i < this.config.sidePoints.length; i++) {
    let pointOption = {
      cx: this.config.sidePoints[i][0],
      cy: this.config.sidePoints[i][1],
      color: "#ffa500",
      r: 4,
      index: parseFloat(i),
      editPoint: true,
    };
    let point = new Points(pointOption, this);
    this.children.push(point);
  }
};
/**
 * make the specific control points darggable based on the given index
 * @param {Number} index the index of the point in the this.control_points
 */
Container.prototype.pointDraggable = function (index) {
  var point_drag = SVG(this.points.nodelist[index]);
  var points = this.points;
  point_drag.draggable();
  // var init_y;
  point_drag.on(`beforedrag.${this.id}point${index}`, () => {
    init_y = points.control_points[0].y;
  });
  point_drag.on(`dragmove.${this.id}point${index}`, (e) => {
    let { handler, box } = e.detail;
    e.preventDefault();
    var x;
    var y;
    switch (index) {
      case 1:
        x = points.control_points[index].x - points.r;
        y = box.cy > points.center_point.y ? points.center_point.y : box.cy;
        break;
      case 5:
        x = points.control_points[index].x - points.r;
        y = box.cy < points.center_point.y ? points.center_point.y : box.cy;
        break;
      case 3:
      case 7:
        x = box.cx;
        y = points.control_points[index].y - points.r;
        break;
      case 0:
      case 2:
        x = box.cx;
        y = box.cy > points.center_point.y ? points.center_point.y : box.cy;
        break;
      case 4:
      case 6:
        x = box.cx;
        y = box.cy < points.center_point.y ? points.center_point.y : box.cy;
        break;
    }
    handler.move(x, y);
    this.updatePointArray(index, x, y);
    this.rotater.updateCor();
    this.producePointsAttribute();
    var width = Math.abs(
      points.control_points[0].x - points.control_points[2].x
    );
    var height = Math.abs(
      points.control_points[0].y - points.control_points[6].y
    );
    this.parent.updateShape(width, height);
    this.parent.adjustTextPosition();
    this.parent.adjustText();
  });
  point_drag.on("dragend.updateCoor", () => {
    this.rotater.option.x = points.control_points[1].x;
    this.rotater.option.y = points.control_points[1].y - 20;
  });
};

/**
 * update the coordinate of each points on polygon
 * @param {Number} sourceIndex the index of the point which trigger this method
 * @param {Number} new_x new X coordinate of the trigger point
 * @param {Number} new_y new Y coordinate of the trigger point
 */
Container.prototype.updatePointArray = function (sourceIndex, new_x, new_y) {
  var rectSVG = SVG(this.config.rectNode);
  var pointsArray = rectSVG.array();
  let border_change = {
    x: [
      [6, 7, 0],
      [2, 3, 4],
    ],
    y: [
      [0, 1, 2],
      [4, 5, 6],
    ],
  };
  let target = {};
  for (let key in border_change) {
    if (!Object.hasOwnProperty.call(border_change, key)) {
      continue;
    }
    var orderArray = [0, 0]; //index 0 取new, index 1 取 2*this.center_point.x/y - new;
    for (let j in border_change[key]) {
      var index = border_change[key][j].indexOf(sourceIndex);
      if (index != -1) {
        orderArray[j] = 1;
        orderArray[1 - j] = -1;
      }
    }
    target[key] = orderArray;
  }
  for (let key in target) {
    if (!Object.hasOwnProperty.call(target, key)) {
      continue;
    }
    var value = key == "x" ? new_x : new_y;
    var ArrayIndex = key == "x" ? 0 : 1;
    var index;
    for (let j in target[key]) {
      if (!Object.hasOwnProperty.call(target[key], j)) {
        continue;
      }
      if (target[key][j]) {
        for (let i in border_change[key][j]) {
          if (!Object.hasOwnProperty.call(border_change[key][j], i)) {
            continue;
          }
          var newCorArray = [null, null];
          index = border_change[key][j][i];
          switch (target[key][j]) {
            case 1:
              this.config.sidePoints[index][ArrayIndex] = value;
              pointsArray[index][ArrayIndex] = value;
              newCorArray[ArrayIndex] = value;
              break;
            case -1:
              let result = 2 * this.config.centerPoint[ArrayIndex] - value;
              this.config.sidePoints[index][ArrayIndex] = result;
              pointsArray[index][ArrayIndex] = result;
              newCorArray[ArrayIndex] = result;
              break;
          }
          this.children[index].updateCoordinate(
            index == sourceIndex,
            this.config.sidePoints[index][0],
            this.config.sidePoints[index][1]
          );
        }
      }
    }
  }
  this.option.width =
    Math.round(
      Math.abs(this.config.sidePoints[0][0] - this.config.sidePoints[2][0]) *
        100
    ) / 100;
  this.option.height =
    Math.round(
      Math.abs(this.config.sidePoints[0][1] - this.config.sidePoints[6][1]) *
        100
    ) / 100;
};

/**
 * create a rotater object
 */
Container.prototype.createRotater = function () {
  let sidePoints = this.config.sidePoints;
  let rotaterOption = {
    x: sidePoints[1][0],
    y: sidePoints[1][1] - 20,
    rotate: this.parent.option.style.rotate,
  };
  var rotater = new Rotater(rotaterOption, this);
  // rotater.draggable();
  return rotater;
};

/**
 * produce the points attribute of container rectangle
 */
Container.prototype.producePointsAttribute = function () {
  let points = "";
  for (let i = 0; i < this.config.sidePoints.length; i++) {
    let x = this.config.sidePoints[i][0];
    let y = this.config.sidePoints[i][1];
    points += `${x} ${y} `;
  }
  return points;
};
/**
 * reshape the container and marker based on the drag point
 * @param {Number} index
 * @param {Number} x
 * @param {Number} y
 */
Container.prototype.reshape = function (index, x, y) {
  x = Math.round(x * 100) / 100;
  y = Math.round(y * 100) / 100;
  this.updatePointArray(index, x, y);
  let points = this.producePointsAttribute();
  this.config.rectNode.setAttribute("points", points);
  this.parent.reshape(this.option.width, this.option.height);
  if (this.children[8]) {
    this.children[8].updateCoordinate(
      this.config.sidePoints[1][0],
      this.config.sidePoints[1][1] - 20
    ); //有无rotater
  }
};

Container.prototype.updateNodePosition = function () {
  let rect = SVG(this.config.rectNode);
  let array = rect.array();
  for (let i = 0; i < array.length; i++) {
    this.config.sidePoints[i][0] = array[i][0];
    this.config.sidePoints[i][1] = array[i][1];
  }
  this.config.centerPoint[0] =
    (this.config.sidePoints[0][0] + this.config.sidePoints[2][0]) / 2;
  this.config.centerPoint[1] =
    (this.config.sidePoints[0][1] + this.config.sidePoints[6][1]) / 2;
};

Container.prototype.updateTransform = function (angel) {
  this.count++;
  let nbox = SVG(this.parent.config.ref).rbox(SVG(this.nodes));
  SVG(this.nodes).rotate(angel, nbox.cx, nbox.cy);
  // this.nodes.setAttribute("transform",`rotate(${rotate},${nbox.cx},${nbox.cy})`);
};

export { Container };
