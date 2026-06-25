import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";
import IbsUtils from "./IbsUtils.js";
import CommonElements from "./CommonElements.js";
import { Domain } from "./Domain.js";

/**
 * Create a basic nucleotide chart, and you can add other components to it
 * @param {{id: String, displayName: String, position: {start: {site: Number, display: String}, end: {site: Number, display: String} },
 * coordinate: {horizontal: {start: String, end: String, isLocked: Boolean,}, vertical: {start: String, isLocked: Boolean}, },
 * style: {align: String, height: Number, fontSize: Number, color: String, gradient:String, texture: {type: String, color: String}},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 * @param {IbsCharts} parent
 */
function Nucleotide(option, parent) {
  CommonElements.ProteinOrNucleotide.call(this, option, parent, "nucleotide");
}

Nucleotide.prototype = Object.create(
  CommonElements.ProteinOrNucleotide.prototype
);
Nucleotide.prototype.constructor = Nucleotide;

Nucleotide.prototype.createDomain = function (option) {
  let domain = new Domain(option, this);
  this.childrenCounts["domain"] += 1;
  this.children.push(domain);
  this.domain.push(domain);
  // console.log(this.domain);
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
 *
 * @param {{id: String, position: {start: {site: Number, display: String}, end: {site: Number, display: String}}, name: String,
 * nameTextStyle: {fontFamily: String, fontSize: Number, fontWeight: String, color: String, rotate: Number, location: String},
 * style: {color: String, texture: {type: String, color: String}, gradient: String},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 */
Nucleotide.prototype.createRibosome = function (option) {
  option = option || {};
  if (!option.id || "" === option.id.trim()) {
    option.id = IbsUtils.getRandomId("Ribosome");
  } else {
    option.id = option.id.trim();
  }
  let randomStartSite = IbsUtils.getRandomNumber(
    this.option.position.start.site,
    this.option.position.end.site - 6
  );
  const defaultOption = {
    id: "",
    position: {
      start: {
        site: 0,
        display: "hide",
      },
      end: {
        site: 0,
        display: "hide",
      },
    },
    displayName: "",
    nameTextStyle: {
      fontFamily: "Arial",
      fontSize: 14,
      color: "#333333",
      rotate: 0,
      location: "hide",
    },
    style: {
      color: "#003f7f",
      texture: {
        type: "",
        color: "#333333",
      },
      gradient: "centerToAllAround",
    },
    borderStyle: {
      color: "#000000",
      size: 0,
      isDash: false,
    },
  };
  option = IbsUtils.supplementToStandardObj(defaultOption, option);
  if (0 === option.position.start.site) {
    if (0 === option.position.end.site) {
      option.position.start.site = randomStartSite;
      option.position.end.site = randomStartSite + 6;
    } else {
      option.position.start.site = parseInt(option.position.end.site) - 6;
    }
  } else {
    if (0 === option.position.end.site) {
      option.position.end.site = parseInt(option.position.start.site) + 6;
    }
  }
  let ribosome = new Ribosome(option, this);
  this.children.push(ribosome);
  this.nodes.appendChild(ribosome.nodes);
  ribosome.draggable(this.config.editable);
  this.parent.config.isUndo = false;
  let operation = {
    target: ribosome,
    cmd: "create",
    args: JSON.stringify(ribosome.option),
    nodesIndex: ribosome.config.nodesIndex,
  };
  this.parent.parent.config.undoStack.push(operation);
  return ribosome;
};

/**
 *
 * @param {{id: String, position: {start: {site: Number, display: String}, end: {site: Number, display: String}}, name: String,
 * nameTextStyle: {fontFamily: String, fontSize: Number, fontWeight: String, color: String, rotate: Number, location: String},
 * style: {color: String, texture: {type: String, color: String}, gradient: String},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 * @param {Nucleotide} parent
 */
function Ribosome(option, parent) {
  this.parent = parent;
  let standard = this.standardizeConfigAndOption(option);
  this.option = standard.option;
  this.config = standard.config;
  this.nodes = this.createNodes();
  let label = this.updateLabel();
  this.nodes.appendChild(label);
}

Ribosome.prototype.createNodes = function () {
  let pStartY = this.parent.option.coordinate.vertical.start;
  let ribosome = IbsUtils.createSvgElement("g", {
    id: this.config.groupId,
  });
  let gradientUrl = this.parent.parent.getGradientUrl(
    this.option.style.gradient,
    this.option.style.color
  );
  let startY = pStartY + this.parent.option.style.height / 2;
  this.config.rx1 = 0.75 * this.config.width;
  // let ry1 = Math.max(rx1 / 2, this.parent.option.style.height);
  this.config.ry1 = Math.max(0.9 * this.parent.option.style.height, 7);
  this.config.rx2 = this.config.width;
  this.config.ry2 = 1.3 * this.config.ry1;
  let d = `M${this.config.startX},${startY} a${this.config.rx1},${this.config.ry1} 0 1,0 ${this.config.width},0 a${this.config.rx2},${this.config.ry2} 0 1,0 -${this.config.width},0`;
  let gradientPathAttr = {
    id: this.config.groupId + "_box",
    d: d,
    style: `fill:url(#${gradientUrl});stroke-width:${this.option.borderStyle.size};stroke:${this.option.borderStyle.color}`,
  };
  if (this.option.borderStyle.isDash) {
    gradientPathAttr["stroke-dasharray"] = "6 4";
  }
  ribosome.appendChild(IbsUtils.createSvgElement("path", gradientPathAttr));
  let textureUrl = this.parent.parent.getTextureUrl(
    this.option.style.texture.type,
    this.option.style.texture.color
  );
  let fillTexture = textureUrl ? `url(#${textureUrl})` : "none";
  let texturePathAttr = {
    id: this.config.groupId + "_textureBox",
    d: d,
    style: `fill:${fillTexture};stroke-width:0;`,
  };
  ribosome.appendChild(IbsUtils.createSvgElement("path", texturePathAttr));
  // let label = this.updateLabel();
  // ribosome.appendChild(label);
  this.createBorder(ribosome);
  return ribosome;
};

/**
 *
 * @param {{id: String, position: {start: {site: Number, display: String}, end: {site: Number, display: String}}, name: String,
 * nameTextStyle: {fontFamily: String, fontSize: Number, fontWeight: String, color: String, rotate: Number, location: String},
 * style: {color: String, texture: {type: String, color: String}, gradient: String},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 */
Ribosome.prototype.update = function (option) {
  let previousOptionJson = JSON.stringify(this.option);
  option = IbsUtils.supplementToStandardObj(this.option, option);
  let standard = this.standardizeConfigAndOption(option);
  this.option = standard.option;
  this.config = standard.config;
  console.log(standard);
  let updateNodes = this.createNodes();
  this.parent.nodes.replaceChild(updateNodes, this.nodes);
  this.nodes = updateNodes;
  this.nodes.appendChild(this.updateLabel()); // Updating the label before causes an exception to be displayed while dragging
  this.draggable(this.parent.config.editable);
  this.parent.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "update",
    args: previousOptionJson,
    nodesIndex: this.config.nodesIndex,
  };
  this.parent.parent.config.undoStack.push(operation);
  return this;
};

Ribosome.prototype.draggable = function (isDraggable = true) {
  if (isDraggable) {
    let ribosome = SVG("#" + this.config.groupId);
    ribosome.draggable();
    ribosome.on("dragstart.namespace", (e) => {
      this.parent.parent.config.isUndo = false;
      let operation = {
        target: this,
        cmd: "drag",
        args: JSON.stringify(this.option),
        nodesIndex: this.config.nodesIndex,
      };
      this.parent.parent.config.undoStack.push(operation);
    });
    ribosome.on("dragmove.namespace", (e) => {
      e.preventDefault();
      const { handler, box } = e.detail;
      let ribsomeBox = SVG("#" + this.config.groupId + "_box");
      let nucleotideRect = SVG("#" + this.parent.config.groupId + "_rect");
      let { x, y } = box;
      y = ribosome.y();
      let isMin = false;
      let isMax = false;
      let minX =
        nucleotideRect.x() -
        (ribsomeBox.x() - ribosome.x()) -
        0.5 * this.config.rx2;
      let maxX = minX + nucleotideRect.width() - this.config.width;
      if (minX > x) {
        x = minX;
        isMax = false;
        isMin = true;
      } else if (maxX < x) {
        x = maxX;
        isMax = true;
        isMin = false;
      } else {
        isMax = false;
        isMin = false;
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
        let pStartX = this.parent.option.coordinate.horizontal.start;
        let rStartX = ribsomeBox.x() + this.config.rx2 / 2;
        this.option.position.start.site =
          this.parent.option.position.start.site +
          parseInt((rStartX - pStartX) / this.parent.config.lengthRatio);
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
    SVG("#" + this.config.groupId).draggable(false);
  }
};

/**
 * update the label with the display: "top" or "bottom", and hide with the display: "hide"
 */
Ribosome.prototype.updateLabel = function () {
  let labelGroupId = this.config.groupId + "_label";
  let labelGroup = document.getElementById(labelGroupId);
  if (!labelGroup) {
    labelGroup = IbsUtils.createSvgElement("g", {
      id: labelGroupId,
    });
  }
  let labelOption = IbsUtils.getDomainLabelOption(this, this.parent);
  let difLength = this.parent.option.style.height / 2;
  // set start label
  let startLabelOption = labelOption.startLabelOption;
  startLabelOption.type = "polyline";
  startLabelOption.y =
    "bottom" === this.option.position.start.display
      ? startLabelOption.y - difLength
      : startLabelOption.y + difLength;
  startLabelOption.lineLength += difLength;
  IbsUtils.updateLabel(startLabelOption, labelGroup);
  // set end label
  let endLabelOption = labelOption.endLabelOption;
  endLabelOption.type = "polyline";
  endLabelOption.y =
    "bottom" === this.option.position.end.display
      ? endLabelOption.y - difLength
      : endLabelOption.y + difLength;
  endLabelOption.lineLength += difLength;
  IbsUtils.updateLabel(endLabelOption, labelGroup);
  return labelGroup;
};

/**
 *
 * @param {{id: String, position: {start: {site: Number, display: String}, end: {site: Number, display: String}}, name: String,
 * nameTextStyle: {fontFamily: String, fontSize: Number, fontWeight: String, color: String, rotate: Number, location: String},
 * style: {color: String, texture: {type: String, color: String}, gradient: String, height: Number},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 */
Ribosome.prototype.standardizeConfigAndOption = function (option) {
  let a = parseInt(option.position.start.site);
  let b = parseInt(option.position.end.site);
  option.position.start.site = Math.max(
    Math.min(a, b),
    this.parent.option.position.start.site
  );
  option.position.end.site = Math.min(
    Math.max(a, b),
    this.parent.option.position.end.site
  );
  option.style.height = parseInt(option.style.height);
  option.style.height = Math.max(
    option.style.height,
    this.parent.option.style.height
  );
  let pStartX = this.parent.option.coordinate.horizontal.start;
  let startX =
    pStartX +
    (option.position.start.site - this.parent.option.position.start.site) *
      this.parent.config.lengthRatio;
  let length = option.position.end.site - option.position.start.site;
  let nodesIndex = this.parent.children.indexOf(this);
  nodesIndex = -1 === nodesIndex ? this.parent.children.length : nodesIndex;
  let config = {
    groupId: this.parent.parent.config.id + "_ribosome_" + option.id,
    length: length,
    startX: startX,
    width: length * this.parent.config.lengthRatio,
    nodesIndex: nodesIndex,
  };
  return {
    option: option,
    config: config,
  };
};

/**
 *
 * @param {SVGElement} ribosome
 */
Ribosome.prototype.createBorder = function (ribosome) {
  if (this.parent.config.editable) {
    let pStartY = this.parent.option.coordinate.vertical.start;
    let x1 = this.config.startX - 0.5 * this.config.rx2;
    let x2 = this.config.startX + 1.5 * this.config.rx2;
    let y1 = pStartY - 1.2 * this.config.ry2;
    let y2 = pStartY + 2.2 * this.config.ry1;
    let borderCoordinates = [
      [x1, y1],
      [x2, y1],
      [x1, y2],
      [x2, y2],
    ];
    let border = IbsUtils.createBorder(borderCoordinates);
    IbsUtils.showBorderByMousedown(ribosome, border);
    ribosome.appendChild(border);
  }
};

Ribosome.prototype.delete = function () {
  IbsUtils.deleteObj(this);
  this.parent.parent.config.isUndo = false;
  let operation = {
    target: this,
    cmd: "delete",
    args: JSON.stringify(this.option),
    nodesIndex: this.config.nodesIndex,
  };
  this.parent.parent.config.undoStack.push(operation);
  return true;
};

export default {
  Nucleotide,
};
