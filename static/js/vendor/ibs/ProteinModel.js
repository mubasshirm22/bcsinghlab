import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";
import IbsUtils from "./IbsUtils.js";
import CommonElements from "./CommonElements.js";

/**
 * Create a basic protein chart, and you can add other components to it
 * @param {{id: String, displayName:String, position: {start: {site: Number, display: String}, end: {site: Number, display: String} },
 * coordinate: {horizontal: {start: String, end: String, isLocked: Boolean,}, vertical: {start: String, isLocked: Boolean}, },
 * style: {align: String, height: Number, fontSize: Number, color: String, gradient:String, texture: {type: String, color: String}},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 * @param {IbsCharts} parent
 */
function Protein(option, parent) {
  CommonElements.ProteinOrNucleotide.call(this, option, parent, "protein");
}

Protein.prototype = Object.create(CommonElements.ProteinOrNucleotide.prototype);
Protein.prototype.constructor = Protein;

/**
 * Create a domain based on one protein
 * @param {{id: String, displayName: String, position: {start: {site: *, display: String}, end: {site: *, display: String}},
 * name: String, nameTextStyle: {fontFamily: String, fontSize: Number, fontWeight: String, color: String, rotate: Number, location: String},
 * style: {color: String, texture: {type: String, color: String}, gradient: String}, borderStyle: {color: String, size: Number, isDash: Boolean}}} domainOption
 */
Protein.prototype.createDomain = function (domainOption) {
  let domain = new Domain(domainOption, this);
  this.nodes.appendChild(domain.nodes);
  domain.draggable(this.config.editable);
  this.children.push(domain);
  this.childrenCounts["domain"] += 1;
  this.parent.config.isUndo = false;
  let operation = {
    target: domain,
    cmd: "create",
    args: JSON.stringify(domain.option),
    nodesIndex: domain.nodesIndex,
  };
  this.parent.config.undoStack.push(operation);
  return domain;
};

/**
 * Create a domain based on one protein
 * @param {{id: String, displayName: String, position: {start: {site: Number, display: String}, end: {site: Number, display: String}}, name: String,
 * nameTextStyle: {fontFamily: String, fontSize: Number, fontWeight: String, color: String, rotate: Number, location: String},
 * style: {color: String, texture: {type: String, color: String}, gradient: String}, borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 * @param {Protein} parent
 */
function Domain(option, parent) {
  this.parent = parent;
  this.nodesIndex = parent.children.length;
  let defaultOption = parent.getDefaultOption("domain");
  option = IbsUtils.supplementToStandardObj(defaultOption, option);
  let standard = this.standardizeConfigAndOption(option);
  this.option = standard.option;
  this.config = standard.config;
  this.nodes = this.createNodes();
}

Domain.prototype.standardizeConfigAndOption = function (
  option,
  updateId = true
) {
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
  option.borderStyle.size = Math.max(0, option.borderStyle.size);
  option.nameTextStyle.fontSize = Math.max(0, option.nameTextStyle.fontSize);
  let PStartX = this.parent.option.coordinate.horizontal.start;
  let DStartX =
    PStartX +
    (option.position.start.site - PStartSite) * this.parent.config.lengthRatio;
  let DWidth =
    (option.position.end.site - option.position.start.site) *
    this.parent.config.lengthRatio;
  let config = {
    groupId: this.parent.config.groupId + "_" + option.id,
    startX: DStartX,
    width: DWidth,
    endX: DStartX + DWidth,
    length: option.position.end.site - option.position.start.site,
    type: "domain",
  };
  return { option: option, config: config };
};

Domain.prototype.createNodes = function () {
  let PY = this.parent.option.coordinate.vertical.start;
  let PHeight = this.parent.option.style.height;
  let domainGroupAttr = {
    type: "domain",
    id: this.config.groupId,
    style: `stroke-width:${this.option.borderStyle.size};`,
  };
  let domainGroup = IbsUtils.createSvgElement("g", domainGroupAttr);
  let gradientUrl = this.parent.parent.getGradientUrl(
    this.option.style.gradient,
    this.option.style.color
  );
  let domainRectAttr = {
    id: this.config.groupId + "_rect",
    x: this.config.startX,
    y: PY,
    // rx: 3,
    width: this.config.width,
    height: PHeight,
    style: `fill:url(#${gradientUrl});stroke:${this.option.borderStyle.color};`,
  };
  if (this.option.borderStyle.isDash) {
    domainRectAttr["stroke-dasharray"] = "6 4";
  }
  let domainRect = IbsUtils.createSvgElement("rect", domainRectAttr);
  domainGroup.appendChild(domainRect);
  this.config.refNodes = domainRect;
  let textureUrl = this.parent.parent.getTextureUrl(
    this.option.style.texture.type,
    this.option.style.texture.color
  );
  let fillTexture = textureUrl ? `url(#${textureUrl})` : "none";
  let textureRectAttr = {
    id: this.config.groupId + "_textureRect",
    x: this.config.startX,
    y: PY,
    rx: 3,
    width: this.config.width,
    height: PHeight,
    style: `fill:${fillTexture};stroke-width:0;`,
  };
  domainGroup.appendChild(IbsUtils.createSvgElement("rect", textureRectAttr));
  if ("" !== this.option.name) {
    let nameTextStyle = this.option.nameTextStyle;
    let halfHeight = this.parent.option.style.height / 2;
    let y = PY + halfHeight;
    let visibility = "visible";
    switch (nameTextStyle.location) {
      case "top":
        y = y - this.parent.option.style.height;
        break;
      case "bottom":
        y = y + this.parent.option.style.height;
        break;
      case "hide":
        visibility = "hidden";
        break;
      default:
        // center
        break;
    }
    let nameStartX = this.config.startX + this.config.width / 2;
    let nameAttr = {
      id: this.config.groupId + "_name",
      style: `font-family:${nameTextStyle.fontFamily};font-style:${nameTextStyle.fontStyle};font-size:${nameTextStyle.fontSize}px;font-weight:${nameTextStyle.fontWeight};fill:${nameTextStyle.color}`,
      "text-anchor": "middle",
      x: nameStartX,
      y: y,
      dy: "0.3em",
      visibility: visibility,
      transform: `rotate(${this.option.nameTextStyle.rotate} ${nameStartX} ${y})`,
    };
    // if (0 == this.option.nameTextStyle.rotate) {
    //     nameAttr["text-anchor"] = "middle";
    // }
    let name = IbsUtils.createSvgElement("text", nameAttr);
    name.textContent = this.option.name;
    domainGroup.appendChild(name);
  }
  domainGroup.appendChild(this.updateLabel());
  this.createBorder(domainGroup);
  this.parent.parent.config.selected = this;
  if (this.parent.config.editable) {
    domainGroup.addEventListener("mousedown", (e) => {
      this.parent.parent.config.selected = this;
      this.parent.parent.config.isHideAll = false;
      this.parent.parent.config.selectedMolecular = this.parent;
    });
  }
  return domainGroup;
};

Domain.prototype.createBorder = function (parentNode) {
  let PY = this.parent.option.coordinate.vertical.start;
  let PHeight = this.parent.option.style.height;
  if (this.parent.config.editable) {
    let coordinates = [
      [this.config.startX, PY],
      [this.config.endX, PY],
      [this.config.startX, PY + PHeight],
      [this.config.endX, PY + PHeight],
    ];
    let border = IbsUtils.createBorder(coordinates);
    IbsUtils.showBorderByMousedown(parentNode, border);
    IbsUtils.showBorderAlone(border);
    parentNode.appendChild(border);
  }
};

/**
 * Update domain with new option
 * @param {{id: String, displayName: String, position: {start: {site: *, display: String}, end: {site: *, display: String}}, name: String,
 * nameTextStyle: {fontFamily: String, fontSize: Number, fontWeight: String, color: String, rotate: Number, location: String},
 * style: {color: String, texture: {type: String, color: String}, gradient: String}, borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 */
Domain.prototype.update = function (option) {
  let previousOptionJson = JSON.stringify(this.option);
  option = IbsUtils.supplementToStandardObj(this.option, option);
  let standard = this.standardizeConfigAndOption(option, false);
  this.option = standard.option;
  this.config = standard.config;
  let updateNodes = this.createNodes();
  this.parent.nodes.replaceChild(updateNodes, this.nodes);
  this.nodes = updateNodes;
  this.draggable(this.parent.config.editable);
  this.parent.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "update",
    args: previousOptionJson,
    nodesIndex: this.nodesIndex,
  };
  this.parent.parent.config.undoStack.push(operation);
  return this;
};

/**
 * Set the domain draggable or not
 * @param {Boolean} isDraggable
 */
Domain.prototype.draggable = function (isDraggable = true) {
  if (isDraggable) {
    // let domainGroup = SVG("#" + this.config.groupId);
    let domainGroup = SVG(this.nodes);
    domainGroup.draggable();
    domainGroup.on(`beforedrag.${this.config.groupId}`, (e) => {
      this.parent.parent.config.isUndo = false;
      let operation = {
        target: this,
        cmd: "drag",
        args: JSON.stringify(this.option),
        nodesIndex: this.nodesIndex,
      };
      this.parent.parent.config.undoStack.push(operation);
    });
    domainGroup.on(`dragmove.${this.config.groupId}`, (e) => {
      e.preventDefault();
      let domainRect = SVG(this.config.refNodes);
      let proteinRect = SVG(this.parent.config.refNodes);
      let isMin = false;
      let isMax = false;
      let domainMinX = proteinRect.x() - (domainRect.x() - domainGroup.x());
      let domainMaxX = domainMinX + proteinRect.width() - domainRect.width();
      const { handler, box } = e.detail;
      let { x, y } = box;
      y = domainGroup.y();
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
      handler.move(x, y);

      if (isMin) {
        this.option.position.start.site =
          this.parent.option.position.start.site;
        this.option.position.end.site =
          this.option.position.start.site + this.config.length;
      } else if (isMax) {
        this.option.position.end.site = this.parent.option.position.end.site;
        this.option.position.start.site =
          this.option.position.end.site - this.config.length;
      } else {
        // let DStartX = Number(domainRect.attr().x);
        let DStartX = domainRect.x();
        let PStartX = this.parent.option.coordinate.horizontal.start;
        this.option.position.start.site =
          this.parent.option.position.start.site +
          Math.round((DStartX - PStartX) / this.parent.config.lengthRatio);
        this.option.position.end.site =
          this.option.position.start.site + this.config.length;
      }
      let startLabel = document.getElementById(
        this.config.groupId + "_start_label"
      );
      if (startLabel) {
        startLabel.textContent = this.option.position.start.site;
      }
      let endLabel = document.getElementById(
        this.config.groupId + "_end_label"
      );
      if (endLabel) {
        endLabel.textContent = this.option.position.end.site;
      }
    });
  } else {
    SVG(this.nodes).draggable(false);
  }
};

/**
 * update the label with the display: "top" or "bottom", and hide with the display: "hide"
 */
Domain.prototype.updateLabel = function () {
  let labelGroupId = this.config.groupId + "_label";
  let labelGroup = document.getElementById(labelGroupId);
  if (!labelGroup) {
    labelGroup = IbsUtils.createSvgElement("g", { id: labelGroupId });
  }
  let labelOption = IbsUtils.getDomainLabelOption(this, this.parent);
  let startLabelOption = labelOption.startLabelOption;
  // startLabelOption["y"] = "bottom" === this.option.position.start.display ? startLabelOption["y"] + 2 : startLabelOption["y"] - 2;
  IbsUtils.updateLabel(startLabelOption, labelGroup);
  let endLabelOption = labelOption.endLabelOption;
  // endLabelOption["y"] = "bottom" === this.option.position.end.display ? endLabelOption["y"] + 2 : endLabelOption["y"] - 2;
  IbsUtils.updateLabel(endLabelOption, labelGroup);
  return labelGroup;
};

Domain.prototype.delete = function () {
  IbsUtils.deleteObj(this);
  this.parent.childrenCounts["domain"] -= 1;
  this.parent.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "delete",
    args: JSON.stringify(this.option),
    nodesIndex: this.nodesIndex,
  };
  this.parent.parent.config.undoStack.push(operation);
  return true;
};

export default { Protein };
