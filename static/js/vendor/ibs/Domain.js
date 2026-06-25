import { Text } from "./Text.js";
import IbsUtils from "./IbsUtils.js";
import { methodsTable } from "./methodsTable.js";
import { Container } from "./Container.js";
import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";

function Domain(domainOption, parent) {
  this.parent = parent;
  this.children = []; //0是container, 1是text
  let { option, config } = this.standarizeOptionAndConfig(domainOption);
  this.option = option;
  this.config = config;
  this.nodesIndex = parent.children.length;
  // this.count = 0;
  this.nodes = this.createNodes();
  this.parent.nodes.appendChild(this.nodes);
  if (this.option.shape != "cylinder") {
    this.children[0].children[8].updateCloneNode();
  }
  this.draggable(this.parent.config.editable);

  this.children[0].updateTransform(this.option.style.rotate);
}

Domain.prototype.standarizeOptionAndConfig = function(
  domainOption,
  updateId = true
) {
  let defaultOption = this.parent.getDefaultOption("domain");
  let option = IbsUtils.supplementToStandardObj(defaultOption, domainOption);
  let config = {};
  if (updateId) {
    option.id = IbsUtils.getRandomId("Domain");
  }
  option.displayName = option.displayName ? option.displayName : option.id;
  option.position.start.site = Math.max(
    1,
    parseInt(option.position.start.site)
  );
  option.position.end.site = Math.max(2, parseInt(option.position.end.site));
  if (option.position.start.site > option.position.end.site) {
    let temp = option.position.end.site;
    option.position.end.site = option.position.start.site;
    option.position.start.site = temp;
  }
  let length = option.position.end.site - option.position.start.site;
  let PStartSite = this.parent.option.position.start.site;
  let PEndSite = this.parent.option.position.end.site;
  if (option.position.start.site < PStartSite) {
    option.position.start.site = PStartSite;
    option.position.end.site = Math.min(PEndSite, PStartSite + length);
  } else if (option.position.end.site > PEndSite) {
    option.position.end.site = PEndSite;
    option.position.start.site = Math.max(PStartSite, PEndSite - length);
  }
  option.height = Math.max(
    this.parent.option.style.height,
    parseInt(option.height)
  );
  option.borderStyle.size = Math.max(0, option.borderStyle.size);
  option.nameTextStyle.fontSize = Math.max(0, option.nameTextStyle.fontSize);

  let PStartX = this.parent.option.coordinate.horizontal.start; //蛋白坐标
  let DStartX =
    PStartX +
    (option.position.start.site - PStartSite) * this.parent.config.lengthRatio; //domain坐标
  let DWidth =
    Math.round(
      (option.position.end.site - option.position.start.site) *
        this.parent.config.lengthRatio *
        100
    ) / 100; //domain宽度
  option.width = DWidth;
  // option.height = option.style.height;
  config = {
    groupId: this.parent.config.groupId + "_" + option.id,
    startX: DStartX, //domain起点坐标
    width: DWidth, //宽度
    endX: DStartX + DWidth, //domain终点坐标
    length: option.position.end.site - option.position.start.site,
    // nodesIndex: this.parent.children.length,
    type: "domain",
    nodePosition: {
      cx: DStartX + DWidth * 0.5,
      cy:
        this.parent.option.coordinate.vertical.start +
        this.parent.option.style.height / 2
    },
    xMin: this.parent.option.coordinate.horizontal.start,
    xMax: this.parent.option.coordinate.horizontal.end
  };
  if (
    option.shape == "circle" ||
    option.shape == "roundRect" ||
    option.shape == "rect"
  ) {
    var { sidePoints, centerPoint } = IbsUtils.calRectEightPointsWithRotate(
      config.nodePosition.cx,
      config.nodePosition.cy,
      option.width,
      option.height,
      option.style.rotate
    );
  } else {
    var { sidePoints, centerPoint } = IbsUtils.calRectEightPoints(
      config.nodePosition.cx,
      config.nodePosition.cy,
      option.width,
      option.height
    );
  }
  config.sidePoints = sidePoints;
  config.centerPoint = centerPoint;
  // config.xMin = this.parent.option.coordinate.horizontal.start + config.width * 0.5;
  // config.xMax = this.parent.option.coordinate.horizontal.end - config.width * 0.5;
  if (option.shape == "cylinder") {
    option.style.rotate = 0;
  }
  return {
    option,
    config
  };
};

Domain.prototype.createNodes = function() {
  // let PY = this.parent.option.coordinate.vertical.start;
  // let PHeight = this.parent.option.style.height;
  let domainGroupAttr = {
    // type: "domain",
    id: this.config.groupId
  };
  let domainGroup = IbsUtils.createSvgElement("g", domainGroupAttr);

  let ref = this.createReferencePoint();
  ref.setAttribute("visibility", "hidden");

  let content = this.createContent();
  this.config.contentNode = content;

  let container = this.createContainer();
  this.children.push(container);

  if (this.option.style.texture.type) {
    var texture = this.createTexture();
    this.config.textureNode = texture;
  }

  if (this.option.name) {
    var text = this.createText();
    this.children.push(text);
  }

  if (this.parent.config.editable) {
    IbsUtils.showBorderByMousedown(domainGroup, container.nodes);
    IbsUtils.showBorderAlone(container.nodes);
  } else {
    IbsUtils.hideBorder(container.nodes);
  }
  domainGroup.appendChild(content);
  if (container) {
    domainGroup.appendChild(container.nodes);
  }
  if (this.option.shape != "cylinder") {
    domainGroup.appendChild(container.children[8].config.cloneNode);
  }
  if (this.config.textureNode) {
    domainGroup.appendChild(this.config.textureNode);
  }
  if (this.children[1]) {
    domainGroup.appendChild(text.nodes);
  }

  domainGroup.appendChild(ref);
  domainGroup.appendChild(this.updateLabel());
  this.parent.parent.config.selected = this;

  domainGroup.addEventListener("mousedown", e => {
    this.parent.parent.config.selected = this;
    this.parent.parent.config.isHideAll = false;
    this.parent.parent.config.selectedMolecular = this.parent;
  });

  return domainGroup;
};

Domain.prototype.createContent = function() {
  let contentOption = methodsTable[this.option.shape].call(
    this,
    this.config.sidePoints
  );
  if (contentOption.points) {
    let pointString = [];
    for (let i = 0; i < contentOption.points.length; i++) {
      let coor = contentOption.points[i];
      let res;
      if (
        this.option.shape == "roundRect" ||
        this.option.shape == "cylinder" ||
        this.option.shape == "circle" ||
        this.option.shape == "rect"
      ) {
        res = contentOption.points[i];
      } else {
        res = IbsUtils.calculatePointAfterRotate(
          this.config.nodePosition.cx,
          this.config.nodePosition.cy,
          coor[0],
          coor[1],
          this.option.style.rotate
        );
      }
      pointString.push(res.join(" "));
    }
    contentOption.points = pointString.join(" ");
  }
  contentOption.id = this.option.id + "_shape";
  let gradientUrl = this.parent.parent.getGradientUrl(
    this.option.style.gradient,
    this.option.style.color
  );
  contentOption.style = `fill:url(#${gradientUrl});stroke:${this.option.borderStyle.color};stroke-width:${this.option.borderStyle.size};`;
  if (this.option.borderStyle.isDash) {
    contentOption["stroke-dasharray"] = "6 4";
  }
  let content = IbsUtils.createShape(this.option.shape, contentOption);
  return content;
};

Domain.prototype.createContainer = function() {
  let containerOption = {
    position: {
      cx: this.config.nodePosition.cx,
      cy: this.config.nodePosition.cy
    },
    width: this.option.width,
    height: this.option.height,
    id: this.option.id + "Container"
  };
  let container = new Container(containerOption, this);
  return container;
};

Domain.prototype.createText = function() {
  let textOption = this.standarizeTextOption(this.option.nameTextStyle);
  textOption.style.borderColor = textOption.style.color;
  textOption.editable = false;
  let text = new Text(textOption, this);
  text.textContent = this.option.name;
  if (this.option.nameTextStyle.location == "hide") {
    text.nodes.setAttribute("visibility", "hidden");
  }
  return text;
};

Domain.prototype.standarizeTextOption = function(option) {
  let textOption = {};
  let position = {};
  var contentBox = SVG(this.config.contentNode).bbox();
  position.x = contentBox.cx;
  switch (option.location.toLowerCase()) {
    case "center":
      position.y = contentBox.cy;
      break;
    case "top":
      position.y = contentBox.cy - contentBox.height * 0.5 - 10;
      break;
    case "bottom":
      position.y = contentBox.cy + contentBox.height * 0.5 + 20;
      break;
    case "hide":
      position.y = contentBox.cy;
      break;
  }
  textOption.content = this.option.name;
  textOption.style = {};
  for (let key in option) {
    if (key == "location") {
      continue;
    }
    textOption.style[key] = option[key];
  }
  textOption.position = position;
  textOption.style.alignmentBaseline = "central";
  return textOption;
};

Domain.prototype.createTexture = function() {
  var texture = this.config.contentNode.cloneNode();
  var textureUrl = this.parent.parent.getTextureUrl(
    this.option.style.texture.type,
    this.option.style.texture.color
  );
  let fillTexture = textureUrl ? `url(#${textureUrl})` : "none";
  texture.style = `stroke-width:0;fill:${fillTexture};`;
  texture.setAttribute("id", `${this.option.id}Texture`);
  return texture;
};

Domain.prototype.createReferencePoint = function() {
  let pointAttr = {
    cx: this.config.nodePosition.cx,
    cy: this.config.nodePosition.cy,
    r: 1,
    fill: "#000000"
  };
  let point = IbsUtils.createSvgElement("circle", pointAttr);
  this.config.ref = point;
  return point;
};

Domain.prototype.reshape = function(width, height, rotate) {
  //更新属性
  if (width) {
    this.option.width = width;
    this.config.width = width;
  }
  // height && (this.option.height = height, this.option.style.height = height);
  if (height) {
    this.option.height = height;
  }
  if (rotate) {
    this.option.style.rotate = rotate;
  }
  let ref = SVG(this.config.ref).bbox();
  // this.calCoorAfterTransform();
  if (
    this.option.shape == "roundRect" ||
    this.option.shape == "circle" ||
    this.option.shape == "rect"
  ) {
    let res = IbsUtils.calRectEightPointsWithRotate(
      this.config.nodePosition.cx,
      this.config.nodePosition.cy,
      this.option.width,
      this.option.height,
      this.option.style.rotate
    );
    this.config.sidePoints = res.sidePoints;
  } else {
    let res = IbsUtils.calRectEightPoints(
      this.config.nodePosition.cx,
      this.config.nodePosition.cy,
      this.option.width,
      this.option.height
    );
    this.config.sidePoints = res.sidePoints;
  }
  //更新图形
  let contentOption = methodsTable[this.option.shape].call(
    this,
    this.config.sidePoints
  );

  if (contentOption.points) {
    let rotatePoints = [];
    if (
      this.option.shape == "roundRect" ||
      this.option.shape == "cylinder" ||
      this.option.shape == "circle" ||
      this.option.shape == "rect"
    ) {
      rotatePoints = contentOption.points;
    } else {
      for (let i = 0; i < contentOption.points.length; i++) {
        let coor = contentOption.points[i];
        let res = IbsUtils.calculatePointAfterRotate(
          this.config.nodePosition.cx,
          this.config.nodePosition.cy,
          coor[0],
          coor[1],
          this.option.style.rotate
        );
        rotatePoints.push(res);
      }
    }
    SVG(this.config.contentNode).plot(rotatePoints);
    SVG(this.config.textureNode).plot(rotatePoints);
  }
  if (!rotate) {
    this.config.length = Math.round(
      this.option.width / this.parent.config.lengthRatio
    );
    let DStartX = ref.cx - this.option.width / 2;
    let PStartX = this.parent.option.coordinate.horizontal.start;
    this.option.position.start.site =
      this.parent.option.position.start.site +
      parseInt((DStartX - PStartX) / this.parent.config.lengthRatio);
    this.option.position.end.site =
      this.option.position.start.site + this.config.length;
    this.config.startX = DStartX;
    this.updateLabel();
  }
  this.adjustTextPosition();
};

Domain.prototype.calCoorAfterTransform = function() {
  let { sidePoints, centerPoint } = IbsUtils.calRectEightPointsWithRotate(
    this.config.nodePosition.cx,
    this.config.nodePosition.cy,
    this.option.width,
    this.option.height,
    this.option.style.rotate
  );
  this.config.sidePoints = sidePoints;
};

Domain.prototype.draggable = function(isDraggable = true) {
  let group = SVG(this.nodes);
  if (isDraggable) {
    let ref = SVG(this.config.ref);
    // let content = SVG(this.config.contentNode);
    let nucleotide = SVG("#" + this.parent.config.groupId + "_rect");
    var gbox;
    var init_cx;
    group.draggable();
    group.on("beforedrag.namespace", e => {
      this.parent.parent.config.isUndo = false;
      let operation = {
        target: this,
        cmd: "drag",
        args: JSON.stringify(this.option),
        nodesIndex: this.nodesIndex
      };
      this.parent.parent.config.undoStack.push(operation);
      gbox = group.bbox();
      init_cx = gbox.cx;
    });
    group.on("dragmove.namespace", e => {
      let { box, handler } = e.detail;
      e.preventDefault();
      let isMin = false;
      let isMax = false;
      // let domainMinX = nucleotide.x() - (content.x() - group.x());
      let domainMinX =
        nucleotide.x() - (ref.cx() - this.option.width / 2 - group.x());
      // let domainMaxX = domainMinX + nucleotide.width() - content.width();
      let domainMaxX = domainMinX + nucleotide.width() - this.option.width;
      let { x, y } = box;
      if (domainMinX > x) {
        x = domainMinX;
        isMin = true;
        isMax = false;
      } else if (domainMaxX < x) {
        x = domainMaxX;
        isMax = true;
        isMin = false;
      } else {
        isMin = false;
        isMax = false;
      }
      group.x(x);

      if (isMin) {
        this.option.position.start.site = this.parent.option.position.start.site;
        this.option.position.end.site =
          this.option.position.start.site + this.config.length;
      } else if (isMax) {
        this.option.position.end.site = this.parent.option.position.end.site;
        this.option.position.start.site =
          this.option.position.end.site - this.config.length;
      } else {
        let DStartX = Number(ref.cx() - this.option.width / 2);
        let PStartX = this.parent.option.coordinate.horizontal.start;
        this.option.position.start.site =
          this.parent.option.position.start.site +
          parseInt((DStartX - PStartX) / this.parent.config.lengthRatio);
        this.option.position.end.site =
          this.option.position.start.site + this.config.length;
      }
      this.config.startX =
        this.config.ref.getAttribute("cx") - this.option.width / 2;
      // console.log(content.width());
      this.updateLabel();
    });
    group.on("dragend.namespance", e => {
      this.children[0].updateNodePosition(ref.bbox().cx, ref.bbox().cy);
      this.updateNodePosition();
    });
  } else {
    group.draggable(false);
  }
};

Domain.prototype.updateTransform = function(x, y, rotate) {
  let angel = rotate - this.option.style.rotate;
  this.option.style.rotate = rotate;
  this.calCoorAfterTransform();
  this.reshape(null, null, rotate);
  this.children[0].updateTransform(angel);
};

Domain.prototype.updateLabel = function() {
  let labelGroupId = this.config.groupId + "_label";
  let labelGroup = document.getElementById(labelGroupId);
  if (!labelGroup) {
    labelGroup = IbsUtils.createSvgElement("g", {
      id: labelGroupId
    });
  }
  let labelOption = IbsUtils.getDomainLabelOption(this, this.parent);
  let startLabelOption = labelOption.startLabelOption;
  startLabelOption.lineLength = (this.option.height * 1) / 2 + 5;
  IbsUtils.updateLabel(startLabelOption, labelGroup);

  let endLabelOption = labelOption.endLabelOption;
  endLabelOption.lineLength = (this.option.height * 1) / 2 + 5;
  IbsUtils.updateLabel(endLabelOption, labelGroup);
  return labelGroup;
};

Domain.prototype.updateNodePosition = function() {
  //获得container矩形边框的中心坐标
  let cx = this.config.ref.getAttribute("cx");
  let cy = this.config.ref.getAttribute("cy");
  this.config.nodePosition.cx = parseFloat(cx);
  this.config.nodePosition.cy = parseFloat(cy);
  this.children[0].updateNodePosition();
};

Domain.prototype.update = function(domainOption) {
  let previousOption = JSON.stringify(this.option);
  this.children = [];
  let defaultOption = this.parent.getDefaultOption("domain");
  domainOption = IbsUtils.supplementToStandardObj(defaultOption, domainOption);
  let originalId = domainOption.id;
  let { option, config } = this.standarizeOptionAndConfig(domainOption, false);
  option.id = originalId;
  this.option = option;
  this.config = config;
  let updateNode = this.createNodes();
  this.parent.nodes.replaceChild(updateNode, this.nodes);
  this.nodes = updateNode;
  this.children[0].updateTransform(this.option.style.rotate);
  if (this.option.shape != "cylinder") {
    this.children[0].children[8].updateCloneNode();
  }
  if (this.parent.config.editable) {
    this.draggable();
  }
  this.parent.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "update",
    args: previousOption,
    nodesIndex: this.nodesIndex
  };
  this.parent.parent.config.undoStack.push(operation);
  return this;
};

Domain.prototype.adjustTextPosition = function() {
  var new_y;
  var contentBox = SVG(this.config.contentNode).bbox();
  switch (this.option.nameTextStyle.location) {
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

Domain.prototype.delete = function() {
  IbsUtils.deleteObj(this);
  // this.parent.childrenCounts["domain"] -= 1;
  this.parent.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "delete",
    args: JSON.stringify(this.option),
    nodesIndex: this.nodesIndex
  };
  this.parent.parent.config.undoStack.push(operation);
  return true;
};

export { Domain };
