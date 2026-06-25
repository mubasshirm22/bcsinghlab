import { Rotater } from "./Rotater.js";
import { Points } from "./Points.js";
import IbsUtils from "./IbsUtils.js";
import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";

function Bracket(bracketOption, parent) {
  this.parent = parent;
  let { option, config } = this.standarizeOption(bracketOption);
  this.children = [];
  this.option = option;
  this.config = config;
  this.nodesIndex = parent.children.length;
  this.nodes = this.createNodes();
  this.parent.nodes.appendChild(this.nodes);
  this.draggable(this.option.editable);
  this.updateTransform();
}

Bracket.prototype.standarizeOption = function(option, updateId = true) {
  let defaultOption = this.parent.getDefaultOption("bracket");
  option = IbsUtils.supplementToStandardObj(defaultOption, option);
  if (updateId) {
    option.id = IbsUtils.getRandomId("Bracket");
  }
  option.displayName = option.displayName ? option.displayName : option.id;
  option.editable = this.parent.config.editable;
  option.style.size = Math.max(1, parseInt(option.style.size));
  option.style.width = Math.max(1, parseInt(option.style.width));
  if (option.shape == "brace") {
    var config = {
      startX: option.coordinate.cx - 10 - option.style.width * 0.5,
      startY: option.coordinate.cy - 10,
      type: "bracket"
    };
    config.d = [
      "M",
      config.startX,
      config.startY,
      "a",
      5,
      5,
      0,
      0,
      0,
      5,
      5,
      "h",
      option.style.width * 0.5,
      "a",
      5,
      5,
      0,
      0,
      1,
      5,
      5,
      "a",
      5,
      5,
      0,
      0,
      1,
      5,
      -5,
      "h",
      option.style.width * 0.5,
      "a",
      5,
      5,
      0,
      0,
      0,
      5,
      -5
    ];
  } else {
    var config = {
      startX: option.coordinate.cx - option.style.width * 0.5,
      startY: option.coordinate.cy - 5,
      type: "bracket"
    };
    config.d = [
      "M",
      config.startX,
      config.startY,
      "v",
      5,
      "h",
      option.style.width,
      "v",
      -5
    ];
  }
  config.nodesTransform = [0, 0, 0, option.coordinate.cx, option.coordinate.cy];
  config.nodeCoor = [option.coordinate.cx, option.coordinate.cy];
  return {
    option,
    config
  };
};

Bracket.prototype.calculateConfig = function(startX) {
  this.config.startX = startX;
  this.config.d[1] = startX;
  this.config.d[12] = this.option.coordinate.cx - startX - 10;
  this.config.d[30] = this.option.coordinate.cx - startX - 10;
  let d = this.config.d.join(" ");
  this.config.bracketNode.setAttribute("d", this.config.d);
};

Bracket.prototype.createNodes = function() {
  let bracketOption = {
    d: this.config.d.join(" "),
    style: `fill:none;stroke:${this.option.style.color};stroke-width:${this.option.style.size}`
  };
  let bracket = IbsUtils.createSvgElement("path", bracketOption);
  let group = IbsUtils.createSvgElement("g", {
    id: `${this.parent.config.id}_${this.option.id}`
  });
  this.config.bracketNode = bracket;
  group.appendChild(bracket);

  let outline = this.createOutline(bracket);
  this.config.outlineNode = outline;
  group.appendChild(outline);
  this.parent.config.selected = this;
  let container = this.createContainer(bracket);
  group.appendChild(container);
  group.appendChild(this.children[2].config.cloneNode);
  if (this.option.editable) {
    IbsUtils.showBorderAlone(this.config.containerNode);
    IbsUtils.showBorderByMousedown(group, this.config.containerNode);
    group.addEventListener("mousedown", () => {
      this.parent.config.selected = this;
    });
  } else {
    IbsUtils.hideBorder(this.config.containerNode);
  }
  return group;
};

Bracket.prototype.createOutline = function(bracket) {
  let outline = bracket.cloneNode();
  outline.style.opacity = 0;
  outline.style.stroke = "black";
  outline.style.strokeWidth = 10;
  return outline;
};

Bracket.prototype.createContainer = function(bracket) {
  let pathArray = SVG(bracket).array();
  // console.log(pathArray);
  let start = pathArray[0].slice(-2);
  let point1Option = {
    cx: start[0],
    cy: start[1],
    r: 4,
    color: "#ffa500",
    index: 0,
    editPoint: true
  };
  var point2Option;
  switch (this.option.shape) {
    case "brace":
      point2Option = {
        cx: pathArray[6][6],
        cy: pathArray[6][7]
      };
      break;
    case "square":
      point2Option = {
        cx: pathArray[2][1],
        cy: pathArray[3][1]
      };
      break;
  }
  point2Option.r = 4;
  point2Option.color = "#ffa500";
  point2Option.index = 1;
  point2Option.editPoint = true;
  let point1 = new Points(point1Option, this);
  let point2 = new Points(point2Option, this);
  this.children.push(point1);
  this.children.push(point2);
  let rotater = this.createRotater();
  this.children.push(rotater);
  let container = IbsUtils.createSvgElement("g", {
    class: "nodeBorder"
  });
  container.appendChild(point1.nodes);
  container.appendChild(point2.nodes);
  container.appendChild(rotater.nodes);
  this.config.containerNode = container;
  return container;
};

Bracket.prototype.createRotater = function() {
  let nbox = SVG(this.config.bracketNode).bbox();
  let rotaterOption = {
    x: parseFloat(nbox.cx) - 6,
    y: parseFloat(nbox.y) + 20,
    rotate: this.option.style.rotate
  };
  let rotater = new Rotater(rotaterOption, this);
  this.children.push(rotater);
  return rotater;
};

Bracket.prototype.draggable = function(isDraggable = true) {
  let group = SVG(this.nodes);
  if (isDraggable) {
    let nbox;
    let dx;
    let dy;
    let targetX;
    let targetY;
    group.draggable();
    group.on("beforedrag.namespace", e => {
      this.parent.config.isUndo = false;
      let operation = {
        target: this,
        cmd: "drag",
        args: JSON.stringify(this.option),
        nodesIndex: this.nodesIndex
      };
      this.parent.config.undoStack.push(operation);
      nbox = SVG(this.nodes).bbox();
      dx = this.option.coordinate.cx - nbox.cx;
      dy = this.option.coordinate.cy - nbox.cy;
    });
    group.on("dragmove.namespace", e => {
      e.preventDefault();
      let { box, handler } = e.detail;
      targetX = Math.round(100 * (parseFloat(box.cx) + dx)) / 100;
      targetY = Math.round(100 * (parseFloat(box.cy) + dy)) / 100;
      this.Move(targetX, targetY);
    });
    group.on("dragend.namespace", e => {
      if (this.children.length) {
        this.children[2].updateCloneNode();
      }
    });
  } else {
    group.draggable(false);
  }
};

Bracket.prototype.Move = function(input_x, input_y) {
  // console.log("drag");
  let translateX =
    typeof input_x == "number"
      ? (translateX = input_x - this.option.coordinate.cx)
      : false;
  let translateY =
    typeof input_y == "number"
      ? (translateY = input_y - this.option.coordinate.cy)
      : false;
  this.option.coordinate.cx = input_x;
  this.option.coordinate.cy = input_y;
  this.updateTransform(translateX, translateY);
};

Bracket.prototype.updateTransform = function(x, y, angel) {
  if (typeof x == "number") {
    this.config.nodesTransform[0] += x;
  }
  if (typeof y == "number") {
    this.config.nodesTransform[1] += y;
  }
  if (typeof angel == "number") {
    this.option.style.rotate = angel;
    this.config.nodesTransform[2] = angel;
  }

  var transformValue = `translate(${this.config.nodesTransform[0]},${this.config.nodesTransform[1]}),rotate(${this.option.style.rotate},${this.config.nodesTransform[3]},${this.config.nodesTransform[4]})`;
  this.config.bracketNode.setAttribute("transform", transformValue);
  this.config.outlineNode.setAttribute("transform", transformValue);
  this.config.containerNode.setAttribute("transform", transformValue);
};

Bracket.prototype.reshape = function(index, startX) {
  if (this.option.shape == "brace") {
    this.config.startX = startX;
    this.config.d[1] = startX;
    this.config.d[12] = this.config.nodeCoor[0] - startX - 10;
    this.config.d[30] = this.config.nodeCoor[0] - startX - 10;
    this.option.style.width = Math.round(this.config.d[12] * 2 * 100) / 100;
  } else {
    this.config.startX = startX;
    this.config.d[1] = startX;
    this.config.d[6] = (this.config.nodeCoor[0] - startX) * 2;
    this.option.style.width =
      Math.round((this.config.nodeCoor[0] - startX) * 200) / 100;
  }
  switch (index) {
    case 0:
      this.children[1].updateCoordinate(
        false,
        2 * this.config.nodeCoor[0] - startX
      );
      break;
    case 1:
      this.children[0].updateCoordinate(false, startX);
      break;
  }
  let d = this.config.d.join(" ");
  this.config.bracketNode.setAttribute("d", d);
  this.config.outlineNode.setAttribute("d", d);
  this.updatePathArray();
};

Bracket.prototype.updatePathArray = function() {
  var dArray = SVG(this.config.bracketNode).array();
  if (this.option.shape == "brace") {
    dArray[0][1] = parseFloat(this.config.d[1]);
    dArray[1][6] = parseFloat(this.config.d[1]) + 5;
    dArray[2][1] =
      parseFloat(this.config.d[1]) + 5 + parseFloat(this.config.d[12]);
    dArray[3][6] =
      parseFloat(this.config.d[1]) + 5 + parseFloat(this.config.d[12]) + 5;
    dArray[4][6] = parseFloat(dArray[3][6]) + 5;
    dArray[5][1] = parseFloat(dArray[4][6]) + parseFloat(this.config.d[30]);
    dArray[6][6] = parseFloat(dArray[5][1]) + 5;
  } else {
  }
};

Bracket.prototype.update = function(bracketOption) {
  let previousOption = JSON.stringify(this.option);
  this.children = [];
  let { option, config } = this.standarizeOption(bracketOption, false);
  this.option = option;
  this.config = config;
  let updateNode = this.createNodes();
  this.parent.nodes.replaceChild(updateNode, this.nodes);
  this.nodes = updateNode;
  this.draggable(this.option.editable);
  this.updateTransform();
  // IbsUtils.showBorderByMousedown(this.nodes, this.config.containerNode);
  // this.nodes.addEventListener("mousedown", () => {
  //   this.parent.config.selected = this;
  // })
  this.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "update",
    args: previousOption,
    nodesIndex: this.nodesIndex
  };
  this.parent.config.undoStack.push(operation);
  return this;
};

Bracket.prototype.delete = function() {
  IbsUtils.deleteObj(this);
  // this.parent.config.elementsCount.bracket -= 1;
  this.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "delete",
    args: JSON.stringify(this.option),
    nodesIndex: this.nodesIndex
  };
  this.parent.config.undoStack.push(operation);
};

export { Bracket };
