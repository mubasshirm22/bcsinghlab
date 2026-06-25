import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";
import IbsUtils from "./IbsUtils.js";
import { Site } from "./Site.js";

/*
 * Common elements of IBS
 */

/**
 * Create a basic protein or nucleotide chart, and you can add other components to it
 * @param {{id: String, displayName:String, position: {start: {site: Number, display: String}, end: {site: Number, display: String} },
 * coordinate: {horizontal: {start: Number, end: Number, isLocked: Boolean,}, vertical: {start: Number, isLocked: Boolean}, },
 * style: {align: String, height: Number, fontSize: Number, color: String, gradient:String, texture: {type: String, color: String}},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} ProteinOrNucleotideOption
 * @param {IbsCharts} parent
 * @param {String} type
 */
function ProteinOrNucleotide(
  ProteinOrNucleotideOption,
  parent,
  type = "protein"
) {
  this.parent = parent;
  this.children = [];
  this.nodesIndex = parent.children.length;
  this.sites = [];
  this.domain = [];
  let standardConfig = this.standardizeConfigAndOption(
    ProteinOrNucleotideOption,
    type
  );
  this.option = standardConfig.option;
  // this.option.coordinate.horizontal.isLocked = !this
  this.config = standardConfig.config;
  this.childrenCounts = {
    domain: 0,
    site: 0,
    cutline: 0,
  };
  this.nodes = this.createNodes();
}

/**
 *
 * @param {{id: String, position: {start: {site: Number, display: String}, end: {site: Number, display: String} },
 * coordinate: {horizontal: {start: Number, end: Number, isLocked: Boolean,}, vertical: {start: Number, isLocked: Boolean}, },
 * style: {align: String, height: Number, fontSize: Number, color: String, gradient:String, texture: {type: String, color: String} },
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 * @param {String} type
 * @param {String} previousOption
 */
ProteinOrNucleotide.prototype.standardizeConfigAndOption = function (
  option,
  type,
  updateId = true
) {
  if (updateId) {
    let upperType = type[0].toUpperCase() + type.slice(1);
    option.id = IbsUtils.getRandomId(upperType);
  }
  option.displayName = option.displayName ? option.displayName : option.id;
  let IbsConfig = this.parent.config;
  // option.coordinate.horizontal.isLocked = !IbsConfig.editable;
  // option.coordinate.vertical.isLocked = !IbsConfig.editable;
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
  option.style.height = Math.max(parseInt(option.style.height), 1);
  option.borderStyle.size = Math.max(0, parseInt(option.borderStyle.size));
  let pOrNLength = option.position.end.site - option.position.start.site;
  if ("custom" === option.style.align) {
    option.coordinate.horizontal.start = parseInt(
      option.coordinate.horizontal.start
    );
    option.coordinate.horizontal.end = parseInt(
      option.coordinate.horizontal.end
    );
  } else if ("rightCustom" === option.style.align) {
    option.coordinate.horizontal.start = parseInt(
      option.coordinate.horizontal.start
    );
    option.coordinate.horizontal.end =
      option.coordinate.horizontal.start +
      parseInt(pOrNLength * IbsConfig.lengthRatio);
  } else if ("selfLocation" === option.style.align) {
    option.coordinate.horizontal.start = parseInt(
      IbsConfig.proteinOrNucleotideStartX +
        (option.position.start.site - 1) * IbsConfig.lengthRatio
    );
    option.coordinate.horizontal.end = parseInt(
      IbsConfig.proteinOrNucleotideStartX +
        (option.position.end.site - 1) * IbsConfig.lengthRatio
    );
  } else if ("justify" === option.style.align) {
    option.coordinate.horizontal.start = IbsConfig.proteinOrNucleotideStartX;
    option.coordinate.horizontal.end = IbsConfig.proteinOrNucleotideEndX;
  } else if ("right" === option.style.align) {
    option.coordinate.horizontal.end = IbsConfig.proteinOrNucleotideEndX;
    option.coordinate.horizontal.start = parseInt(
      option.coordinate.horizontal.end - pOrNLength * IbsConfig.lengthRatio
    );
  } else {
    // the default align is "left"
    option.style.align = "left";
    option.coordinate.horizontal.start = IbsConfig.proteinOrNucleotideStartX;
    option.coordinate.horizontal.end =
      IbsConfig.proteinOrNucleotideStartX +
      parseInt(pOrNLength * IbsConfig.lengthRatio);
  }
  let Width =
    option.coordinate.horizontal.end - option.coordinate.horizontal.start;
  let originalRatio = Width / pOrNLength;
  let lengthRatio =
    originalRatio < 0.1
      ? Number(originalRatio.toPrecision(2))
      : Number(originalRatio.toFixed(2));
  let config = {
    lengthRatio: lengthRatio,
    width: Width,
    length: pOrNLength,
    groupId: IbsConfig.id + "_" + option.id,
    editable: IbsConfig.editable,
    type: type,
    nodeBorder: null,
  };

  return {
    option: option,
    config: config,
  };
};

ProteinOrNucleotide.prototype.createNodes = function () {
  let groupAttr = {
    id: this.config.groupId,
    // type: "Protein",
    style: `font-size:${this.option.style.fontSize}px;stroke-width:${this.option.borderStyle.size};`,
  };
  let group = IbsUtils.createSvgElement("g", groupAttr);
  let StartY = this.option.coordinate.vertical.start;
  let gradientUrl = this.parent.getGradientUrl(
    this.option.style.gradient,
    this.option.style.color
  );
  let width =
    this.option.coordinate.horizontal.end -
    this.option.coordinate.horizontal.start;
  let RectAttr = {
    id: this.config.groupId + "_rect",
    x: this.option.coordinate.horizontal.start,
    y: StartY,
    width: width,
    height: this.option.style.height,
    style: `fill:url(#${gradientUrl});stroke:${this.option.borderStyle.color};`,
  };
  if (this.option.borderStyle.isDash) {
    RectAttr["stroke-dasharray"] = "6 4";
  }
  let Rect = IbsUtils.createSvgElement("rect", RectAttr);
  group.appendChild(Rect);
  this.config.refNodes = Rect;
  let textureUrl = this.parent.getTextureUrl(
    this.option.style.texture.type,
    this.option.style.texture.color
  );
  let fillTexture = textureUrl ? `url(#${textureUrl})` : "none";
  let textureRectAttr = {
    id: this.config.groupId + "_textureRect",
    x: this.option.coordinate.horizontal.start,
    y: StartY,
    width: width,
    height: this.option.style.height,
    style: `fill:${fillTexture};stroke-width:0;`,
  };
  group.appendChild(IbsUtils.createSvgElement("rect", textureRectAttr));
  let label = this.updateLabel();
  group.appendChild(label);
  this.createBorder(group);
  this.parent.config.selected = this;
  this.parent.config.selectedMolecular = this;
  if (this.config.editable) {
    group.addEventListener("mousedown", (e) => {
      this.parent.config.selected = this;
      this.parent.config.isHideAll = false;
      this.parent.config.selectedMolecular = this;
    });
  }
  return group;
};

ProteinOrNucleotide.prototype.createBorder = function (targetElement) {
  if (this.config.editable) {
    let borderCoordinates = [
      [
        this.option.coordinate.horizontal.start,
        this.option.coordinate.vertical.start - 9,
      ],
      [
        this.option.coordinate.horizontal.end,
        this.option.coordinate.vertical.start - 9,
      ],
      [
        this.option.coordinate.horizontal.start,
        this.option.coordinate.vertical.start + this.option.style.height + 9,
      ],
      [
        this.option.coordinate.horizontal.end,
        this.option.coordinate.vertical.start + this.option.style.height + 9,
      ],
    ];
    let border = IbsUtils.createBorder(borderCoordinates);
    this.config.nodeBorder = border;
    IbsUtils.showBorderByMousedown(targetElement, border);
    IbsUtils.showBorderAlone(border);
    targetElement.appendChild(border);
  }
};

/**
 * update the label with the display: "top" or "bottom", and hide with the display: "hide"
 */
ProteinOrNucleotide.prototype.updateLabel = function () {
  let labelGroupId = this.config.groupId + "_label";
  let labelGroup = document.getElementById(labelGroupId);
  if (!labelGroup) {
    labelGroup = IbsUtils.createSvgElement("g", {
      id: labelGroupId,
    });
  }
  // set the start label
  let labelLength = Math.max(6, this.option.style.height / 2);
  let startLabelOption = {
    id: this.config.groupId + "_start",
    type: "line",
    position: "start",
    display: this.option.position.start.display,
    label: this.option.position.start.site,
    x: this.option.coordinate.horizontal.start,
    lineLength: labelLength,
  };

  startLabelOption["y"] =
    "bottom" === this.option.position.start.display
      ? this.option.coordinate.vertical.start + this.option.style.height
      : this.option.coordinate.vertical.start;
  IbsUtils.updateLabel(startLabelOption, labelGroup);
  // set the end label
  let endLabelOption = {
    id: this.config.groupId + "_end",
    type: "line",
    position: "end",
    display: this.option.position.end.display,
    label: this.option.position.end.site,
    x: this.option.coordinate.horizontal.end,
    lineLength: labelLength,
  };
  endLabelOption["y"] =
    "bottom" === this.option.position.end.display
      ? this.option.coordinate.vertical.start + this.option.style.height
      : this.option.coordinate.vertical.start;
  IbsUtils.updateLabel(endLabelOption, labelGroup);
  return labelGroup;
};

/**
 *
 * @param {{horizontal: Boolean, vertical: Boolean}} params
 */
ProteinOrNucleotide.prototype.draggable = function (
  params = {
    horizontal: true,
    vertical: true,
  }
) {
  // let id = this.config.groupId;
  let group = SVG(this.nodes);
  // this.option.coordinate.horizontal.isLocked = true === params.horizontal ? false : true;
  this.option.coordinate.horizontal.isLocked = !params.horizontal;
  // this.option.coordinate.vertical.isLocked = true === params.vertical ? false : true;
  this.option.coordinate.vertical.isLocked = !params.vertical;
  let horizontalLocked = this.option.coordinate.horizontal.isLocked;
  let verticalLocked = this.option.coordinate.vertical.isLocked;

  if (horizontalLocked && verticalLocked) {
    group.draggable(false);
  } else {
    group.draggable();
    group.on(`beforedrag.${this.config.groupId}`, (e) => {
      this.parent.config.isUndo = false;
      let operation = {
        target: this,
        cmd: "drag",
        args: JSON.stringify(this.option),
        nodesIndex: this.nodesIndex,
      };
      this.parent.config.undoStack.push(operation);
    });
    group.on(`dragmove.namespace.${this.config.groupId}`, (e) => {
      const { handler, box } = e.detail;
      e.preventDefault();
      let { x, y } = box;
      x = horizontalLocked ? group.x() : x;
      y = verticalLocked ? group.y() : y;
      handler.move(x, y);
    });
    group.on(`dragend.${this.config.groupId}`, (e) => {
      this.option.style.align = "rightCustom";
      let groupRect = SVG(this.config.refNodes);
      this.option.coordinate.horizontal.start = parseInt(groupRect.x());
      this.option.coordinate.horizontal.end =
        this.option.coordinate.horizontal.start + this.config.width;
      this.option.coordinate.vertical.start = parseInt(groupRect.y());
      for (let i = 0; i < this.sites.length; i++) {
        this.sites[i].updateNodePosition();
        for (let j = 0; j < this.sites[i].config.siteLines.length; j++) {
          let line = this.sites[i].config.siteLines[j];
          var dArray = SVG(line.config.lineNode).array();
          line.updatePathArray(0, dArray[0][1], dArray[0][2]);
          for (let k = 1; k < dArray.length; k++) {
            line.updatePathArray(k, dArray[k][1], dArray[k][2]);
          }
        }
      }
      //更新domain的中心点
      if (this.config.type == "nucleotide") {
        for (let i = 0; i < this.domain.length; i++) {
          this.domain[i].updateNodePosition();
        }
      }
    });
  }
};

ProteinOrNucleotide.prototype.delete = function () {
  IbsUtils.deleteObj(this);
  // this.parent.config.elementsCount["proteinOrNucleotide"] -= 1;
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
 * Update
 * @param {{id: String, displayName: String, position: {start: {site: Number, display: String}, end: {site: Number, display: String} },
 * coordinate: {horizontal: {start: String, end: String, isLocked: Boolean,}, vertical: {start: String, isLocked: Boolean}, },
 * style: {align: String, height: Number, fontSize: Number, color: String, gradient:String, texture: {type: String, color: String}},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 */
ProteinOrNucleotide.prototype.update = function (option) {
  let previousOption = JSON.stringify(this.option);
  option = IbsUtils.supplementToStandardObj(this.option, option);
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
  let pConfig = this.parent.config;
  if (
    1 === pConfig.elementsCount.proteinOrNucleotide &&
    length !== pConfig.scaleLength
  ) {
    pConfig.scaleLength = length;
    let originalRatio =
      (pConfig.proteinOrNucleotideEndX - pConfig.proteinOrNucleotideStartX) /
      length;
    pConfig.lengthRatio =
      originalRatio < 0.1
        ? Number(originalRatio.toPrecision(2))
        : Number(originalRatio.toFixed(2));
    // Math.round(100 * ((pConfig.proteinOrNucleotideEndX - pConfig.proteinOrNucleotideStartX) / length)) / 100;
  }
  let standard = this.standardizeConfigAndOption(
    option,
    this.config.type,
    false
  );
  this.option = standard.option;
  this.config = standard.config;
  let updateProteinOrNucleotide = this.createNodes();
  for (let i = 0; i < this.children.length; i++) {
    // let updateChild = this.children[i].update(this.children[i].option, true)
    let updateChild = this.children[i].update(this.children[i].option);
    updateProteinOrNucleotide.appendChild(updateChild.nodes);
    this.parent.config.undoStack.pop();
  }
  this.parent.nodes.replaceChild(updateProteinOrNucleotide, this.nodes);
  this.nodes = updateProteinOrNucleotide;
  IbsUtils.showBorderAlone(this.config.nodeBorder);
  this.draggable({
    horizontal: !this.option.coordinate.horizontal.isLocked,
    vertical: !this.option.coordinate.vertical.isLocked,
  });
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

ProteinOrNucleotide.prototype.getDefaultOption = function (type = "domain") {
  // let id = IbsUtils.getRandomId(type);
  let defaultOption;
  let pHeight = this.option.style.height;
  let defaultLength = Math.max(1, Math.floor(this.config.length / 20));
  switch (type) {
    case "domain":
      let defaultDomainLength =
        10 >= this.config.length ? 1 : Math.floor(this.config.length / 10);
      let randomDomainStartSite = IbsUtils.getRandomNumber(
        this.option.position.start.site,
        this.option.position.end.site - defaultDomainLength
      );
      defaultOption = {
        id: IbsUtils.getRandomId("Domain"),
        position: {
          isLocked: false,
          start: {
            site: randomDomainStartSite,
            display: "hide",
          },
          end: {
            site: randomDomainStartSite + defaultDomainLength,
            display: "hide",
          },
        },
        displayName: `Domain${this.childrenCounts[type] + 1}`,
        text: "",
        nameTextStyle: {
          fontFamily: "Arial",
          fontStyle: "normal",
          fontSize: 14,
          color: "#333333",
          rotate: 0,
          location: "center",
          fontWeight: "normal",
          rotate: 0,
        },
        style: {
          color: "#0000ff",
          texture: {
            type: "none",
            color: "#4389b4",
          },
          gradient: "centerToTopAndBottom",
        },
        borderStyle: {
          color: "#000000",
          size: 1,
          isDash: false,
        },
      };
      switch (this.config.type) {
        case "nucleotide":
          defaultOption.shape = "roundRect";
          defaultOption.height = 5 * pHeight;
          defaultOption.style.color = "#4389b4";
          // defaultOption.style.height = 5 * pHeight;
          defaultOption.style.rotate = 0;
          break;
        default:
          // protein
          break;
      }
      break;
    case "cutline":
      const randomStartSite = IbsUtils.getRandomNumber(
        this.option.position.start.site,
        this.option.position.end.site - defaultLength
      );
      const defaultHeight =
        "protein" === this.config.type ? 1.5 * pHeight : 3 * pHeight;
      defaultOption = {
        id: IbsUtils.getRandomId("Cutline"),
        displayName: `Cutline${this.childrenCounts[type] + 1}`,
        position: {
          start: {
            site: randomStartSite,
          },
          end: {
            site: randomStartSite + defaultLength,
          },
        },
        style: {
          height: defaultHeight,
          type: "diagonal",
          direction: "forward",
        },
        borderStyle: {
          color: "#000000",
          size: 1,
          isDash: false,
        },
      };
      break;
    case "site":
      let randomSite = IbsUtils.getRandomNumber(
        this.option.position.start.site,
        this.option.position.end.site
      );
      defaultOption = {
        id: IbsUtils.getRandomId("Site"),
        displayName: `Site${this.childrenCounts[type] + 1}`,
        shape: "circle",
        mode: "site",
        location: randomSite,
        coordinate: {
          cy: this.option.coordinate.vertical.start - 50,
          cx:
            Math.round(
              100 *
                ((randomSite - this.option.position.start.site) *
                  this.config.lengthRatio)
            ) /
              100 +
            this.option.coordinate.horizontal.start,
        },
        width: 25,
        height: 25,
        text: {
          content: "",
          position: "top",
          style: {
            fontSize: 14,
            fontFamily: "Arial",
            color: "#333333",
            // location: "up",
            fontStyle: "normal",
            fontWeight: "normal",
            rotate: 0,
          },
        },
        style: {
          color: "#efc326",
          borderColor: "#000000",
          borderSize: 1,
          isDash: false,
          rotate: 0,
          gradient: "centerToAllAround",
          texture: {
            type: "none",
            color: "#4389b4",
          },
        },
        editable: this.config.editable,
      };
      break;
    default:
      break;
  }

  return defaultOption;
};

/**
 * Create a cutline on the protein or nucleotide
 * @param {{id: String, displayName:String, position: {start: {site: Number}, end: {site: Number},} , style: { height: Number, type: String, direction: String}, borderStyle: {color: String, size: Number}}} option
 */
ProteinOrNucleotide.prototype.createCutline = function (option) {
  let defaultOption = this.getDefaultOption("cutline");
  let pHeight = this.option.style.height;
  option = IbsUtils.supplementToStandardObj(defaultOption, option);
  if (pHeight > option.style.height || 5 * pHeight < option.style.height) {
    option.style.height = defaultOption.style.height;
  }
  let tempEndSite = option.position.end.site;
  if (option.position.start.site > tempEndSite) {
    option.position.end.site = option.position.start.site;
    option.position.start.site = tempEndSite;
  }
  option.position.start.site = Math.max(
    this.option.position.start.site,
    option.position.start.site
  );
  option.position.end.site = Math.min(
    this.option.position.end.site,
    option.position.end.site
  );
  let cutline = new CutLine(option, this);
  this.children.push(cutline);
  this.childrenCounts["cutline"] += 1;
  this.nodes.appendChild(cutline.nodes);
  cutline.draggable(this.config.editable);
  this.parent.config.isUndo = false;
  let operation = {
    target: cutline,
    cmd: "create",
    args: JSON.stringify(cutline.option),
    nodesIndex: cutline.config.nodesIndex,
  };
  this.parent.config.undoStack.push(operation);
  return cutline;
};

ProteinOrNucleotide.prototype.createSite = function (option) {
  let site = new Site(option, this);
  // console.log(site);
  this.children.push(site);

  this.sites.push(site);
  this.childrenCounts.site += 1;

  // console.log(site);
  this.parent.config.isUndo = false;
  let operation = {
    target: site,
    cmd: "create",
    args: JSON.stringify(site.option),
    nodesIndex: site.nodesIndex,
  };
  this.parent.config.undoStack.push(operation);
  return site;
};

/**
 * 生成 Cutline
 * @param {{id: String, displayName:String, position: {start: {site: Number}, end: {site: Number},} , style: { height: Number, type: String, direction: String}, borderStyle: {color: String, size: Number}}} option
 * @param {ProteinOrNucleotide} parent
 */
function CutLine(option, parent) {
  this.parent = parent;
  let standard = this.standardizeConfigAndOption(option);
  this.option = standard.option;
  this.config = standard.config;
  this.nodes = this.createNodes();
  this.nodesIndex = parent.children.length;
}

CutLine.prototype.createNodes = function () {
  let group = IbsUtils.createSvgElement("g", {
    id: this.config.groupId,
  });
  // let parentStartX = this.parent.option.coordinate.horizontal.start;
  // let startX = parentStartX + (this.option.position.start.site - this.parent.option.position.start.site) * this.parent.config.lengthRatio;
  let startX = this.config.startX;
  let startY = this.config.startY;
  // let width = (this.option.position.end.site - this.option.position.start.site) * this.parent.config.lengthRatio;
  let singalWidth = 0.8 * this.config.width;
  let intervalWidth = 0.2 * this.config.width;
  let d1 = ""; // the fill path
  let d2 = ""; // the border line path
  switch (this.option.style.type) {
    case "diagonal":
      let l1 = "";
      let l2 = "";
      if ("opposite" === this.option.style.direction) {
        startY += this.option.style.height;
        l1 = `l${singalWidth},-${this.option.style.height}`;
        l2 = `l-${singalWidth},${this.option.style.height}`;
      } else {
        l1 = `l${singalWidth},${this.option.style.height}`;
        l2 = `l-${singalWidth},-${this.option.style.height}`;
      }
      d1 = `M${startX},${startY} ${l1} h${intervalWidth} ${l2}`;
      d2 = `M${startX},${startY} ${l1} m${intervalWidth},0 ${l2}`;
      break;
    case "bezier":
      let x1 = 0.25 * singalWidth;
      let y1 = 0.75 * this.option.style.height;
      let x2 = 0.75 * singalWidth;
      let y2 = 0.25 * this.option.style.height;
      let c1 = "";
      let c2 = "";
      if ("opposite" === this.option.style.direction) {
        startY += this.option.style.height;
        c1 = `c${x1},-${y1} ${x2},-${y2} ${singalWidth},-${this.option.style.height}`;
        c2 = `c-${x1},${y1} -${x2},${y2} -${singalWidth},${this.option.style.height}`;
      } else {
        // forward
        c1 = `c${x1},${y1} ${x2},${y2} ${singalWidth},${this.option.style.height}`;
        c2 = `c-${x1},-${y1} -${x2},-${y2} -${singalWidth},-${this.option.style.height}`;
      }
      d2 = `M${startX},${startY} ${c1} m${intervalWidth},0 ${c2}`;
      d1 = `M${startX},${startY} ${c1} h${intervalWidth} ${c2}`;
      break;
    default:
      // diagonal
      break;
  }
  let fillAttr = {
    d: d1,
    style: `stroke-width:0;fill:white;`,
  };
  let fillPath = IbsUtils.createSvgElement("path", fillAttr);
  let lineAttr = {
    d: d2,
    style: `stroke:${this.option.borderStyle.color};fill:none;stroke-width:${this.option.borderStyle.size};`,
  };
  if (this.option.borderStyle.isDash) {
    lineAttr["stroke-dasharray"] = "3 2";
  }
  let linePath = IbsUtils.createSvgElement("path", lineAttr);

  group.appendChild(fillPath);
  group.appendChild(linePath);
  this.createBorder(group);
  this.parent.parent.config.selected = this;
  if (this.parent.config.editable) {
    group.addEventListener("mousedown", (e) => {
      this.parent.parent.config.selected = this;
      this.parent.parent.config.isHideAll = false;
      this.parent.parent.config.selectedMolecular = this.parent;
    });
  }
  return group;
};

/**
 * standardize thr cutline config and option
 * @param {{id: String, displayName:String, position: {start: {site: Number}, end: {site: Number},} , style: { height: Number, type: String, direction: String}, borderStyle: {color: String, size: Number}}} option
 */
CutLine.prototype.standardizeConfigAndOption = function (
  option,
  updateId = true
) {
  if (updateId) {
    option.id = IbsUtils.getRandomId("Cutline");
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
  option.style.height = Math.max(
    parseInt(option.style.height),
    this.parent.option.style.height
  );
  option.borderStyle.size = Math.max(0, parseInt(option.borderStyle.size));
  let parentStartX = this.parent.option.coordinate.horizontal.start;

  let difHeight = (option.style.height - this.parent.option.style.height) / 2;
  let startY = this.parent.option.coordinate.vertical.start - difHeight;
  let config = {
    groupId: this.parent.config.groupId + "_" + option.id,
    startX:
      parentStartX +
      (option.position.start.site - this.parent.option.position.start.site) *
        this.parent.config.lengthRatio,
    width: length * this.parent.config.lengthRatio,
    length: option.position.end.site - option.position.start.site,
    startY: startY, // Upper left coordinates
    type: "cutline",
  };
  return {
    option: option,
    config: config,
  };
};

CutLine.prototype.createBorder = function (parentNode) {
  if (this.parent.config.editable) {
    let endX = this.config.startX + this.config.width;
    let endY = this.config.startY + this.option.style.height;
    let coordinates = [
      [this.config.startX, this.config.startY],
      [endX, this.config.startY],
      [this.config.startX, endY],
      [endX, endY],
    ];
    let border = IbsUtils.createBorder(coordinates);
    IbsUtils.showBorderByMousedown(parentNode, border);
    IbsUtils.showBorderAlone(border);
    parentNode.appendChild(border);
  }
};

CutLine.prototype.update = function (option) {
  let previousOptionJson = JSON.stringify(this.option);
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

CutLine.prototype.draggable = function (isDraggable = true) {
  let cutline = SVG("#" + this.config.groupId);
  if (isDraggable) {
    cutline.draggable();
    cutline.on("dragstart.namespace", (e) => {
      this.parent.parent.config.isUndo = false;
      let operation = {
        target: this,
        cmd: "drag",
        args: JSON.stringify(this.option),
        nodesIndex: this.nodesIndex,
      };
      this.parent.parent.config.undoStack.push(operation);
    });
    cutline.on("dragmove.namespace", (e) => {
      e.preventDefault();
      let proteinRect = SVG("#" + this.parent.config.groupId + "_rect");
      let minX = proteinRect.x() - 4;
      let maxX = minX + proteinRect.width() - this.config.width;
      let isMin = false;
      let isMax = false;
      const { handler, box } = e.detail;
      let { x, y } = box;
      y = cutline.y();
      if (x < minX) {
        x = minX;
        isMin = true;
        isMax = false;
      } else if (x > maxX) {
        x = maxX;
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
        this.config.startX = this.parent.option.coordinate.horizontal.start;
      } else if (isMax) {
        this.option.position.end.site = this.option.position.end.site;
        this.option.position.start.site =
          this.option.position.end.site - this.config.length;
        this.config.startX = cutline.x() + 4;
      } else {
        let startX = cutline.x() + 4;
        let pStartX = this.parent.option.coordinate.horizontal.start;
        this.config.startX = startX;
        this.option.position.start.site =
          this.parent.option.position.start.site +
          parseInt((startX - pStartX) / this.parent.config.lengthRatio);
        this.option.position.end.site =
          this.option.position.start.site + this.config.length;
      }
    });
  } else {
    cutline.draggable(false);
  }
};

/**
 * Delete this object
 */
CutLine.prototype.delete = function () {
  IbsUtils.deleteObj(this);
  // this.parent.childrenCounts.cutline -= 1;
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

export default {
  ProteinOrNucleotide,
  CutLine,
};
