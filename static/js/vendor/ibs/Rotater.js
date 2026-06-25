import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";
import IbsUtils from "./IbsUtils.js";

/**
 * create a circle node used to rotate the parent charts
 * @param {{x:Number, y:Number, rotate: Number}} option rotaterOption
 * @param {{IbsCharts}} parent
 */
function Rotater(option, parent) {
  this.option = option;
  this.parent = parent; //对应group元素
  this.config = {
    type: "rotater",
  };
  this.nodes = this.createNodes();
  this.cloneRotater();
  this.draggable();
}
/**
 * create a rotater node
 */
Rotater.prototype.createNodes = function () {
  let rotaterAttr = {
    d: `M ${this.option.x} ${this.option.y} a 5 5 0 1 1 5 -5 v 3 M ${
      this.option.x + 2
    } ${this.option.y - 5} l 3 4 l 3 -4`,
    fill: "none",
    stroke: "red",
    "stroke-width": "2",
  };
  var node = IbsUtils.createSvgElement("path", rotaterAttr);
  this.config.d = [
    ["M", this.option.x, this.option.y],
    ["a", 5, 5, 0, 1, 1, 5, -5],
    ["v", 3],
    ["M", this.option.x + 2, this.option.y - 5],
    ["l", 3, 4],
    ["l", 3, -4],
  ];
  return node;
};
/**
 * create a transparent rotater
 */
Rotater.prototype.cloneRotater = function () {
  let cloneNode = this.nodes.cloneNode();
  // let nbox = SVG(this.parent.nodes).bbox();
  let d = this.config.d.slice(0);
  d[0][1] = 0;
  d[0][2] = 0;
  d[3][1] = d[0][1] + 2;
  d[3][2] = d[0][2] - 5;
  let dAttr = [];
  for (let i = 0; i < d.length; i++) {
    dAttr.push(d[i].join(" "));
  }
  cloneNode.setAttribute("d", dAttr.join(" "));
  cloneNode.setAttribute("id", `${this.parent.option.id}CloneRotater`);
  cloneNode.style.fill = "white";
  cloneNode.style.opacity = "0";
  cloneNode.style.visibility = "visible";
  cloneNode.style.stroke = "black";
  cloneNode.style["stroke-width"] = "2";
  this.config.cloneNode = cloneNode;
};
/**
 * make to rotater draggable and rotate the charts
 */
Rotater.prototype.draggable = function () {
  var pointSVG = SVG(this.config.cloneNode);
  pointSVG.draggable();
  var x0;
  var y0;
  var x1;
  var y1;
  var box;
  let parentType = this.parent.config.type;
  pointSVG.on(`beforedrag.${this.parent.option.id}Rotate`, (e) => {
    let IbsRoot;
    let target;
    if ("container" === parentType) {
      // marker, site, domain
      if ("marker" === this.parent.parent.config.type) {
        target = this.parent.parent;
        IbsRoot = target.parent;
      } else {
        // domain, site
        target = this.parent.parent;
        IbsRoot = target.parent.parent;
      }
    } else {
      // text, bracket
      target = this.parent;
      IbsRoot = target.parent;
    }

    IbsRoot.config.isUndo = false;
    let operation = {
      target: target,
      cmd: "rotate",
      args: JSON.stringify(target.option),
      // nodesIndex
    };
    IbsRoot.config.undoStack.push(operation);
  });
  pointSVG.on(`dragstart.${this.parent.option.id}Rotate`, (e) => {
    if (parentType == "text") {
      box = SVG(this.parent.config.textNode).rbox(pointSVG);
    }
    if (parentType == "container") {
      box = SVG(this.parent.config.rectNode).rbox(pointSVG);
    }
    if (parentType == "bracket") {
      box = SVG(this.parent.config.bracketNode).rbox(pointSVG);
    }
    if (parentType == "line") {
      box = {
        cx: this.option.x,
        cy: this.option.y,
      };
    }
    x0 = parseFloat(box.cx);
    y0 = parseFloat(box.cy);
    x1 = x0;
    //根据父节点不同更改参考点
    if (this.parent.config.type == "bracket") {
      y1 = 20000;
    } else {
      y1 = -20000;
    }
  });
  pointSVG.on(`dragmove.${this.parent.id}Rotate`, (e) => {
    let rotate = this.option.rotate;
    e.preventDefault();
    let { box, handler } = e.detail;

    let angel = IbsUtils.Cosines(
      x1,
      y1,
      parseFloat(box.cx),
      parseFloat(box.cy),
      x0,
      y0
    );
    if (isNaN(angel)) {
      // console.log("NaN");
      return false;
    }
    rotate = Math.round(angel);
    if (rotate > 180) {
      rotate -= 360;
    }
    if (rotate < -180) {
      rotate += 360;
    }
    handler.move(box.cx, box.cy);
    if (this.parent.config.type == "bracket") {
      this.option.rotate = -rotate;
    } else {
      this.option.rotate = rotate;
    }
    if (
      this.parent.config.type == "text" ||
      this.parent.config.type == "bracket"
    ) {
      this.parent.updateTransform(null, null, this.option.rotate);
    } else if (this.parent.config.type == "container") {
      this.parent.parent.updateTransform(null, null, this.option.rotate);
    } else {
      //parentType == line
      this.parent.rotate(this.option.rotate);
    }
  });
  pointSVG.on(`dragend.test`, () => {
    this.updateCloneNode();
  });
};
/**
 * update the coordinate of the clone node
 */
Rotater.prototype.updateCloneNode = function () {
  let gbox = SVG(this.nodes).rbox(SVG(this.parent.parent.nodes));
  let d = this.config.d.slice(0);
  d[0][1] = parseFloat(gbox.x) + parseFloat(gbox.width) / 2;
  d[0][2] = parseFloat(gbox.y) + parseFloat(gbox.height);
  d[3][1] = d[0][1] + 2;
  d[3][2] = d[0][2] - 5;
  SVG(this.config.cloneNode).plot(d);
  // if (this.config.cloneNode.transform) {
  //   this.config.cloneNode.setAttribute("transform", "");
  // }
};
/**
 * update the coordinate of the rotater node
 * @param {Number} x
 * @param {Number} y
 */
Rotater.prototype.updateCoordinate = function (x, y) {
  //console.log("kmr")
  this.option.x = x;
  this.option.y = y;
  this.config.d[0][1] = x;
  this.config.d[0][2] = y;
  this.config.d[3][1] = x + 2;
  this.config.d[3][2] = y - 5;
  SVG(this.nodes).plot(this.config.d);
  this.updateCloneNode();
};

export { Rotater };
