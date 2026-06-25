import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";
import { Text } from "./Text.js";
import IbsUtils from "./IbsUtils.js";
import { methodsTable } from "./methodsTable.js";
import { Container } from "./Container.js";
import { Rotater } from "./Rotater.js";

/**
 * create a marker object
 * @param {{shape:String, x:Number, y:Number, style:{fill:String, stroke:String, strokeWidth: Number},
 * text:{content:String, style: {textAnchor: String}}, editable: Boolean} MarkerOption
 */
function Markers(MarkerOption, parent) {
  this.parent = parent;
  this.children = [];
  let { option, config } = this.standarizeOptionAndConfig(MarkerOption);
  this.option = option;
  this.config = config;
  this.nodesIndex = parent.children.length;
  this.nodes = this.createNodes();
  this.updateTransform();
  this.draggable(this.option.editable);
}

/**
 * standarize the input marker option
 * @param {{}} MarkerOption
 */
Markers.prototype.standarizeOptionAndConfig = function (
  MarkerOption,
  updateId = true
) {
  var config = { type: "marker" };
  let defaultMarker = this.parent.getDefaultOption("marker");
  let option = IbsUtils.supplementToStandardObj(defaultMarker, MarkerOption);
  if (updateId) {
    option.id = IbsUtils.getRandomId("Marker");
  }
  option.displayName = option.displayName ? option.displayName : option.id;
  option.editable = this.parent.config.editable;
  option.coordinate.cx = parseInt(option.coordinate.cx);
  option.coordinate.cy = parseInt(option.coordinate.cy);
  option.width = Math.max(1, parseInt(option.width));
  option.height = Math.max(1, parseInt(option.height));
  option.text.style.fontSize = Math.max(0, option.text.style.fontSize);
  option.style.borderSize = Math.max(0, parseInt(option.style.borderSize));
  let { cx, cy } = option.coordinate;
  let { sidePoints, centerPoint } = IbsUtils.calRectEightPoints(
    option.coordinate.cx,
    option.coordinate.cy,
    option.width,
    option.height
  );
  config.sidePoints = sidePoints;
  config.centerPoint = centerPoint;
  config.nodePosition = {
    cx,
    cy,
  };
  config.nodesTransform = [0, 0, 0];
  return {
    option,
    config,
  };
};
/**
 * create nodes based on the marker option
 */
Markers.prototype.createNodes = function () {
  let content = this.createContent();
  let container;
  // let rotater;
  let group = IbsUtils.createSvgElement("g", {
    id: `${this.option.id}Group`,
  });
  this.parent.config.selected = this;
  // if (this.option.editable) {
  container = this.createContainer();
  this.children.push(container);
  group.addEventListener("mousedown", () => {
    this.parent.config.selected = this;
    this.parent.config.isHideAll = false;
  });
  // }
  if (this.option.style.texture.color) {
    var texture = this.createTexture();
    this.config.textureNode = texture;
  }
  if (this.option.text.content) {
    var text = this.createText();
    this.children.push(text);
  }
  /* 
  // for (let i = 0; i < this.config.siteLines.length; i++) {
    console.log("this.config.siteLines");
    console.log(this.config);
  for (let i in this.config.siteLines) {
    group.appendChild(this.config.siteLines[i].nodes);
  } 
 */

  group.appendChild(content);
  if (texture) {
    group.appendChild(texture);
  }
  if (text) {
    group.appendChild(text.nodes);
  }
  if (container) {
    group.appendChild(container.nodes);
  }
  if (this.option.shape != "cylinder") {
    if (container) {
      group.appendChild(container.children[8].config.cloneNode);
    }
  }

  if (this.option.editable) {
    IbsUtils.showBorderByMousedown(group, container.nodes);
    IbsUtils.showBorderAlone(container.nodes);
  } else {
    IbsUtils.hideBorder(container.nodes);
  }
  return group;
};
/**
 * create the content node
 */
Markers.prototype.createContent = function () {
  let contentOption = methodsTable[this.option.shape].call(
    this,
    this.config.sidePoints
  );
  if (contentOption.points) {
    let pointString = [];
    for (let i = 0; i < contentOption.points.length; i++) {
      let onePoint = contentOption.points[i].join(" ");
      pointString.push(onePoint);
    }
    contentOption.points = pointString.join(" ");
  }
  contentOption.id = this.option.id;
  let content;
  switch (this.option.shape) {
    case "triangle": //
    case "rhombus": //
    case "parallelogram": //
    case "trapezium": //
    case "pentagon": //
    case "hexagon": //
    case "octagon": //
    case "arrow": //
    case "dovetailArrow": //
      content = IbsUtils.createSvgElement("polygon", contentOption);
      break;
    case "roundRect":
    case "rect":
      content = IbsUtils.createSvgElement("rect", contentOption);
      break;
    case "circle":
      content = IbsUtils.createSvgElement("ellipse", contentOption);
      break;
    case "none":
      content = IbsUtils.createSvgElement("circle", contentOption);
      break;
    case "cylinder":
      content = IbsUtils.createSvgElement("path", contentOption);
      break;
    default:
      throw "unsupport shape!";
  }
  let styleOption = IbsUtils.standarizeStyleOption(this.option.style);
  if (this.option.shape != "none") {
    if (this.option.style.gradient !== "none") {
      var gradientUrl = this.parent.getGradientUrl(
        this.option.style.gradient,
        this.option.style.color
      );
      styleOption.fill = `url(#${gradientUrl})`;
    }
  } else {
    styleOption.fill = "#114514";
    styleOption.stroke = "none";
  }
  for (let key in styleOption) {
    content.style[key] = styleOption[key];
  }
  this.config.contentNode = content;
  return content;
};
/**
 * create a container object
 */
Markers.prototype.createContainer = function () {
  let containerOption = {
    position: {
      cx: this.option.coordinate.cx,
      cy: this.option.coordinate.cy,
    },
    width: this.option.width,
    height: this.option.height,
    id: this.option.id + "Container",
  };
  let container = new Container(containerOption, this);
  return container;
};
/**
 * create a rotater object
 */
Markers.prototype.createRotater = function () {
  let sidePoints = this.children[0].config.sidePoints;
  let rotaterOption = {
    x: sidePoints[1][0] - 6,
    y: sidePoints[1][1] - 20,
    rotate: this.option.style.rotate,
  };
  var rotater = new Rotater(rotaterOption, this);
  rotater.draggable();
  return rotater;
};
/**
 * create a text object
 */
Markers.prototype.createText = function () {
  let textOption = this.standarizeTextPosition(this.option.text);
  textOption.style.borderColor = textOption.style.color;
  textOption.id = `${this.option.id}Text`;
  textOption.editable = false;
  let text = new Text(textOption, this);
  if (text.option.location == "center") {
    this.adjustTextContent(text);
  }
  if (this.option.text.position == "hide") {
    text.nodes.setAttribute("visibility", "hidden");
  }
  return text;
};
/**
 * clone the content node used to contain the texture.
 */
Markers.prototype.createTexture = function () {
  var texture = this.config.contentNode.cloneNode();
  var textureUrl = this.parent.getTextureUrl(
    this.option.style.texture.type,
    this.option.style.texture.color
  );
  let fillTexture = textureUrl ? `url(#${textureUrl})` : "none";
  texture.style = `stroke-width:0;fill:${fillTexture};`;
  texture.setAttribute("id", `${this.option.id}Texture`);
  return texture;
};
/**
 *
 */
Markers.prototype.createGradient = function () {
  var gradient = this.config.contentNode.cloneNode();
  var gradientUrl = this.parent.getGradientUrl(
    this.option.style.gradient,
    this.option.style.color
  );
  gradient.style = `stroke:${this.option.style.borderColor};fill:url(#${gradientUrl})`;
  gradient.setAttribute("id", `${this.option.id}Gradient`);
  return gradient;
};

/**
 * update the nodes based on the marker option
 * @param {*} MarkerOption
 */
Markers.prototype.update = function (MarkerOption) {
  let previousOption = JSON.stringify(this.option);
  this.children = [];
  let { option, config } = this.standarizeOptionAndConfig(MarkerOption, false);
  this.option = option;
  this.config = config;
  let updateNode = this.createNodes();
  this.parent.nodes.replaceChild(updateNode, this.nodes);
  this.nodes = updateNode;
  this.updateTransform();
  if (this.option.editable) {
    this.draggable();
  }

  this.children[0].nodes.setAttribute("visibility", "hidden");
  // IbsUtils.showBorderByMousedown(this.nodes,this.children[0].nodes);
  this.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "update",
    args: previousOption,
    nodesIndex: this.nodesIndex,
  };
  this.parent.config.undoStack.push(operation);
  return this;
};
/**
 * delete this marker object
 */
Markers.prototype.delete = function () {
  IbsUtils.deleteObj(this);
  // this.parent.config.elementsCount -= 1;
  this.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "delete",
    args: JSON.stringify(this.option),
    nodesIndex: this.nodesIndex,
  };
  this.parent.config.undoStack.push(operation);
  return true;
};
/**
 * reshape the content, including its width, height and type
 * @param {Number} width content width
 * @param {Number} height content height
 * @param {Number} type content type
 */
Markers.prototype.reshape = function (width, height, shape) {
  this.option.width = width;
  this.option.height = height;
  if (shape) {
    this.option.shape = shape;
  }
  let { sidePoints, centerPoint } = IbsUtils.calRectEightPoints(
    this.config.nodePosition.cx,
    this.config.nodePosition.cy,
    width,
    height
  );
  let contentOption = methodsTable[this.option.shape].call(this, sidePoints);
  if (contentOption.points) {
    let pointString = [];
    for (let i = 0; i < contentOption.points.length; i++) {
      let onePoint = contentOption.points[i].join(" ");
      pointString.push(onePoint);
    }
    contentOption.points = pointString.join(" ");
  }
  for (let key in contentOption) {
    this.config.contentNode.setAttribute(key, contentOption[key]);
    if (this.config.textureNode) {
      this.config.textureNode.setAttribute(key, contentOption[key]);
    }
  }
  if (this.option.text.content) {
    if (this.option.text.location != "center") {
      this.adjustTextPosition();
    } else {
      this.adjustTextContent(this.children[1]);
    }
  }
  this.config.sidePoints = sidePoints;
  this.config.centerPoint = centerPoint;
};

/**
 * transform the position attribute of text option to coordiante
 * @param {{}} option IBS option
 */
Markers.prototype.standarizeTextPosition = function (option) {
  let textOption = {};
  // var positionArray = option[key].split(",");
  let position = {};
  var contentBox = SVG(this.config.contentNode).bbox();
  position.x = contentBox.cx;
  switch (option.position.toLowerCase()) {
    case "center":
      position.y = contentBox.cy;
      // console.log(1);
      break;
    case "top":
      position.y = contentBox.cy - contentBox.height * 0.5 - 10;
      // console.log(2);
      break;
    case "bottom":
      position.y = contentBox.cy + contentBox.height * 0.5 + 20;
      // console.log(3);
      break;
    case "hide":
      position.y = contentBox.cy;
      break;
  }
  textOption.content = option.content;
  textOption.style = {};
  for (let key in option.style) {
    textOption.style[key] = option.style[key];
  }
  textOption.position = position;
  textOption.style.alignmentBaseline = "central";
  return textOption;
};
/**
 * adjust the text position while reshape the content
 */
Markers.prototype.adjustTextPosition = function () {
  var new_y;
  var contentBox = SVG(this.config.contentNode).bbox();
  switch (this.option.text.position) {
    case "top":
      new_y = contentBox.cy - contentBox.height * 0.5 - 10;
      this.children[1].config.textNode.setAttribute("y", new_y);
      break;
    case "bottom":
      new_y = contentBox.cy + contentBox.height * 0.5 + 20;
      this.children[1].config.textNode.setAttribute("y", new_y);
      break;
    default:
      break;
  }
};

/**
 * adjust text content while changing the box size, replace with "..."
 */
Markers.prototype.adjustTextContent = function (textObj) {
  if (this.option.text.location == "center") {
    var cloneText = textObj.config.textNode.cloneNode();
    cloneText.textContent = textObj.option.content;
    var textWidth = SVG(cloneText).bbox().width;
    var textLength = textObj.option.content.length;
    //console.log(textWidth / textLength);
    var perChart = Math.ceil(textWidth / textLength);
    var maxNum = Math.floor(this.option.width / perChart);
    if (textLength > maxNum) {
      var adjustText = textObj.option.content.substring(0, maxNum - 1) + "...";
      textObj.updateContent(adjustText);
    } else {
      textObj.updateContent(textObj.option.content);
    }
  }
};
/**
 * make the marker draggable
 */
Markers.prototype.draggable = function (isDraggable = true) {
  var groupSVG = SVG(this.nodes);
  if (isDraggable) {
    groupSVG.draggable();
    var disX;
    var disY;
    groupSVG.on(`beforedrag.${this.id}Group`, (e) => {
      this.parent.config.isUndo = false;
      let operation = {
        target: this,
        cmd: "drag",
        args: JSON.stringify(this.option),
        nodesIndex: this.nodesIndex,
      };
      this.parent.config.undoStack.push(operation);
      var nbox = IbsUtils.calAbsoluteCoordinate(
        this.children[0].config.rectNode,
        this.nodes
      );
      var mbox = SVG(this.nodes).bbox();
      disX = parseFloat(nbox.cx) - parseFloat(mbox.cx);
      disY = parseFloat(nbox.cy) - parseFloat(mbox.cy);
    });
    groupSVG.on(`dragmove.${this.id}Group`, (e) => {
      let { box } = e.detail;
      e.preventDefault();
      let new_x = Math.round((parseFloat(box.cx) + disX) * 100) / 100;
      let new_y = Math.round((parseFloat(box.cy) + disY) * 100) / 100;
      this.Move(new_x, new_y);
    });
    groupSVG.on("dragend.test", (e) => {
      if (this.option.shape != "cylinder") {
        if (this.children[0]) {
          this.children[0].children[8].updateCloneNode();
        }
      }
    });
  } else {
    groupSVG.draggable(false);
  }
};
/**
 *
 * @param {Number} input_x the direction x coordiante
 * @param {Number} input_y the direction y coordiante
 */
Markers.prototype.Move = function (input_x, input_y, selfMove = true) {
  let nbox = IbsUtils.calAbsoluteCoordinate(
    this.children[0].config.rectNode,
    this.nodes
  );
  var new_x = parseFloat(input_x) - parseFloat(nbox.cx);
  var new_y = parseFloat(input_y) - parseFloat(nbox.cy);
  this.updateTransform(new_x, new_y, null);
  this.option.coordinate.cx = input_x;
  this.option.coordinate.cy = input_y;
};
/**
 * update the transform attribute of the nodes
 * @param {Number} x
 * @param {Number} y
 * @param {Number} angel
 */
Markers.prototype.updateTransform = function (x, y, angel) {
  if (typeof x == "number") {
    this.config.nodesTransform[0] += x;
  }
  if (typeof y == "number") {
    this.config.nodesTransform[1] += y;
  }
  if (typeof angel == "number") {
    this.option.style.rotate = angel;
  }
  var nbox = SVG(this.config.contentNode).bbox();
  var transformValue = `translate(${this.config.nodesTransform[0]},${this.config.nodesTransform[1]}),rotate(${this.option.style.rotate},${nbox.cx},${nbox.cy})`;
  this.config.contentNode.setAttribute("transform", transformValue);
  if (this.children[0]) {
    this.children[0].nodes.setAttribute("transform", transformValue); //有无container
  }
  // this.children[1] && this.children[1].nodes.setAttribute("transform", transformValue); //有无rotater
  if (this.option.text.content) {
    if (this.config.type == "site" && this.option.shape == "none") {
      this.children[1].Move(
        this.option.coordinate.cx,
        this.option.coordinate.cy
      );
    } else {
      this.children[1].updateTransform(x, y);
    }
  }
  if (this.option.style.texture.color) {
    this.config.textureNode.setAttribute("transform", transformValue);
  }
};

Markers.prototype.updateStartEndSite = function () {
  let centerX =
    (this.option.coordinate.cx -
      this.parent.option.coordinate.horizontal.start) /
    this.parent.config.lengthRatio;
  let disX = this.option.width / this.parent.config.lengthRatio;
  this.config.start.site = parseInt(centerX - disX / 2);
  this.config.end.site = parseInt(centerX + disX / 2);
  this.config.center = (this.config.start.site + this.config.end.site) / 2;
  this.option.position.start.site = parseInt(centerX - disX / 2);
  this.option.position.end.site = parseInt(centerX + disX / 2);
  this.config.startX = this.option.coordinate.cx - this.option.width / 2;
  this.config.width = this.option.width;
};
Markers.prototype.updateLabel = function () {
  let { startLabelOption, endLabelOption } = IbsUtils.getDomainLabelOption(
    this,
    this.parent
  );
  let startGroup = IbsUtils.createSvgElement("g", {
    id: this.config.groupId + "_start_label",
  });
  let endGroup = IbsUtils.createSvgElement("g", {
    id: this.config.groupId + "_end_label",
  });
  IbsUtils.updateLabel(startLabelOption, startGroup);
  IbsUtils.updateLabel(endLabelOption, endGroup);
  this.config.label.start = startGroup;
  this.config.label.end = endGroup;
  return {
    startGroup,
    endGroup,
  };
};

export default Markers;
