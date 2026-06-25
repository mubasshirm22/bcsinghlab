import IbsUtils from "./IbsUtils.js";
import ProteinModel from "./ProteinModel.js";
import NucleotideModel from "./NucleotideModel.js";
import { Text } from "./Text.js";
import Marker from "./Marker.js";
import Line from "./Line.js";
import { Bracket } from "./Bracket.js";

/**
 * Based on the prepared DOM, initialize the IBS Object
 * @param {HTMLElement} dom
 * @param {Boolean} editable
 */
function init(dom, editable = true) {
  let config = createIbsConfig(dom, editable);
  return new IbsCharts(dom, config);
}

/**
 * Create a IbsCharts instance
 * @param {HTMLElement} parent
 * @param {{id: String, canvasWidth: Number, canvasHeight: Number, canvasMinWidth: Number, canvasMinHeight: Number, viewBox: String,
 * showCanvasSize: Boolean, showCanvasBorder: Boolean, canvasBorder: *, proteinOrNucleotideStartX: Number, proteinOrNucleotideEndX: Number, defaultLength: Number,
 * scaleLength: Number, lengthRatio: Number, proteinDefaultHeight: Number, nucleotideDefaultHeight: Number, editable: Boolean,
 * curProteinOrNucleotideY: Number, gradientIdList: String[], textureIdList: String[], zoomRatio: Number, zoomRange: [Number, Number], isUndo: Boolean,
 * undoStack: {target:Object, cmd: String, args: String, nodesIndex:Number}[], redoStack: {target:Object, cmd: String, args: String, nodesIndex:Number}[]},
 * selected: *, selectedMolecular: *, type: String} config
 */
function IbsCharts(parent, config) {
  this.config = config;
  this.parent = parent;
  this.children = [];
  let canvasNote = createCanvasNote(config);
  this.parent.appendChild(canvasNote);
  this.noteNode = canvasNote;
  let svgRootAttr = {
    xmlns: "http://www.w3.org/2000/svg",
    // "xmlns:xlink":"http://www.w3.org/1999/xlink", // 当拖动文字时会自动加上该属性，这里也加上会导致下载的SVG因属性重复而报错无法打开
    id: config.id,
    width: config.canvasWidth,
    height: config.canvasHeight,
    viewBox: config.viewBox
  };
  let svgChart = IbsUtils.createSvgElement("svg", svgRootAttr);
  this.defsNode = IbsUtils.createSvgElement("defs", {
    id: config.id + "_defs"
  });
  svgChart.appendChild(this.defsNode);
  parent.appendChild(svgChart);
  this.nodes = svgChart;
  this.showCanvasBorder(this.config.showCanvasBorder);
  this.draggable(this.config.editable);

  this.zoomByMouseWheel(this.config.editable);

  this.hideSelectedBorder();
}

/**
 * show or remove the canvas border
 * @param {Boolean} isShow
 */
IbsCharts.prototype.showCanvasBorder = function(isShow = true) {
  this.config.showCanvasBorder = isShow;
  let canvasBorderId = this.config.id + "_canvasBorder";
  // let border = document.getElementById(canvasBorderId);
  if (isShow) {
    if (this.config.canvasBorder) {
      this.config.canvasBorder.setAttribute("visibility", "visible");
    } else {
      let canvasRectAttr = {
        id: canvasBorderId,
        x: 0,
        y: 0,
        width: this.config.canvasWidth,
        height: this.config.canvasHeight,
        fill: "none",
        stroke: "gray",
        "stroke-width": 1,
        "stroke-dasharray": "3 2",
        visibility: "visible"
      };
      this.config.canvasBorder = IbsUtils.createSvgElement(
        "rect",
        canvasRectAttr
      );
      IbsUtils.insertAfter(this.config.canvasBorder, this.defsNode);
    }
  } else {
    if (this.config.canvasBorder) {
      this.config.canvasBorder.setAttribute("visibility", "hidden");
    }
  }
};

/**
 *
 * @param {Boolean} draggable
 */
IbsCharts.prototype.draggable = function(draggable = true) {
  if (draggable) {
    this.nodes.onmousedown = e => {
      e.preventDefault();
      let preViewBoxList = this.config.viewBox.split(" ");
      let onDragStartX = e.clientX;
      let onDragStartY = e.clientY;
      this.nodes.onmousemove = ev => {
        let moveX = parseInt(preViewBoxList[0]) - (ev.clientX - onDragStartX);
        let moveY = parseInt(preViewBoxList[1]) - (ev.clientY - onDragStartY);
        let curViewBox = `${moveX} ${moveY} ${preViewBoxList[2]} ${preViewBoxList[3]}`;
        this.nodes.setAttribute("viewBox", curViewBox);
        this.config.viewBox = curViewBox;
      };
    };
    this.nodes.onmouseup = e => {
      e.preventDefault();
      this.nodes.onmousemove = null;
    };
    this.nodes.onmouseleave = e => {
      e.preventDefault();
      this.nodes.onmousemove = null;
    };
  } else {
    this.nodes.onmousemove = null;
    this.nodes.onmousedown = null;
    this.nodes.onmouseup = null;
    this.nodes.onmouseleave = null;
  }
};

/**
 * Adjust the canvas to specified size and center it
 * @param {[Number, Number]} param
 */
IbsCharts.prototype.adjustCanvasSize = function() {
  // console.log(this);
  // let preViewBox = this.config.viewBox;
  let preRatio = this.config.zoomRatio;
  // let preStartX = this.config.proteinOrNucleotideStartX;
  // let preEndX = this.config.proteinOrNucleotideEndX;
  this.centerCanvas();
  let chartsWidth = this.config.maxX - this.config.minX;
  // console.log(chartsWidth);
  let chartsHeight = this.config.maxY - this.config.minY;
  // console.log(chartsHeight);
  let ratioX = (0.8 * this.config.canvasWidth) / chartsWidth;
  let ratioY = (0.8 * this.config.canvasHeight) / chartsHeight;
  // console.log(ratioX + " | " + ratioY);
  let ratio = Math.round(Math.min(ratioX, ratioY) * preRatio * 10) / 10;
  // console.log(ratio);
  if (preRatio > ratio) {
    this.zoom(ratio);
  }
};

/**
 * Reset the viewBox to its initial state
 */
IbsCharts.prototype.resetViewBox = function() {
  this.config.viewBox = `0 0 ${this.config.canvasWidth} ${this.config.canvasHeight}`;
  this.nodes.setAttribute("viewBox", this.config.viewBox);
};

IbsCharts.prototype.centerCanvas = function() {
  if (0 < this.children.length) {
    let pOrNMinXList = [];
    let pOrNMaxXList = [];
    let minXList = [];
    let maxXList = [];
    let yList = [];
    for (let i = 0; i < this.children.length; i++) {
      switch (this.children[i].config.type) {
        case "protein":
        case "nucleotide":
          pOrNMinXList.push(
            this.children[i].option.coordinate.horizontal.start
          );
          pOrNMaxXList.push(this.children[i].option.coordinate.horizontal.end);
          yList.push(this.children[i].option.coordinate.vertical.start);
          break;
        case "marker":
          minXList.push(this.children[i].option.coordinate.cx);
          maxXList.push(this.children[i].option.coordinate.cx);
          yList.push(this.children[i].option.coordinate.cy);
          break;
        case "line":
          let coordinates = this.children[i].option.position.split(/\s+/);
          let startPoint = coordinates[0].split(",");
          let endPoiont = coordinates[coordinates.length - 1].split(",");
          let lineMinX = Math.min(
            parseInt(startPoint[0]),
            parseInt(endPoiont[0])
          );
          let lineMaxX = Math.max(
            parseInt(startPoint[0]),
            parseInt(endPoiont[0])
          );
          minXList.push(lineMinX);
          maxXList.push(lineMaxX);
          yList.push(parseInt(startPoint[1]), parseInt(endPoiont[1]));
          break;
        case "text":
          minXList.push(this.children[i].option.position.x);
          maxXList.push(this.children[i].option.position.x);
          yList.push(this.children[i].option.position.y);
          break;
        default:
          break;
      }
    }
    let pOrNMaxX;
    let pOrNMinX;
    if (0 < pOrNMinXList.length) {
      pOrNMaxX = Math.max(...pOrNMaxXList);
      pOrNMinX = Math.min(...pOrNMinXList);
      maxXList.push(pOrNMaxX);
      minXList.push(pOrNMinX);
    }
    this.config.minX = Math.min(...minXList);
    this.config.maxX = Math.max(...maxXList);
    this.config.minY = Math.min(...yList);
    this.config.maxY = Math.max(...yList);
    // console.log("minX | maxX | minY | maxY");
    // console.log(`${minX} | ${maxX} | ${minY} | ${maxY}`);
    let ratio = this.config.zoomRatio;
    let offsetX = Math.round(
      (this.config.maxX + this.config.minX - this.config.canvasWidth * ratio) /
        2
    );
    let offsetY = Math.round(
      (this.config.maxY + this.config.minY - this.config.canvasHeight * ratio) /
        2
    );
    let realOffsetX = Math.round(
      offsetX + (this.config.canvasWidth * ratio - this.config.canvasWidth) / 2
    );
    let realOffsetY = Math.round(
      offsetY +
        (this.config.canvasHeight * ratio - this.config.canvasHeight) / 2
    );
    // this.zoom(1);
    let pOrNOffsetX = pOrNMinX
      ? (pOrNMaxX + pOrNMinX - this.config.canvasWidth * ratio) / 2
      : offsetX;

    // if (pOrNMinX) {
    //     pOrNOffsetX = (pOrNMaxX + pOrNMinX - this.config.canvasWidth * ratio) / 2;
    // }
    let viewBox = `${realOffsetX} ${realOffsetY} ${this.config.canvasWidth} ${this.config.canvasHeight}`;
    this.config.viewBox = viewBox;
    this.nodes.setAttribute("viewBox", viewBox);
    let startX = Math.round(0.1 * this.config.canvasWidth * ratio);
    this.config.proteinOrNucleotideStartX = startX + pOrNOffsetX;
    this.config.proteinOrNucleotideEndX =
      this.config.proteinOrNucleotideStartX +
      Math.round(this.config.lengthRatio * this.config.scaleLength);
    // this.zoom(ratio);
    if (this.config.canvasBorder) {
      this.config.canvasBorder.setAttribute("x", realOffsetX);
      this.config.canvasBorder.setAttribute("y", realOffsetY);
      this.config.canvasBorder.setAttribute("width", this.config.canvasWidth);
      this.config.canvasBorder.setAttribute("height", this.config.canvasHeight);
    }
  }
};

/**
 * zoom into specified point, and set the zoom config with the option
 * @param {Number} ratio
 * @param {{zoomMin: Number, zoomMax: Number}} option
 */
IbsCharts.prototype.zoom = function(
  ratio,
  option = { zoomMin: 0.2, zoomMax: 5 }
) {
  let zoomMin =
    parseFloat(option.zoomMin) <= 0 ? 0.2 : parseFloat(option.zoomMin);
  this.config.zoomRange = [zoomMin, parseFloat(option.zoomMax)];
  let preRatio = this.config.zoomRatio;
  this.config.zoomRatio = parseFloat(ratio);
  let changeRatio = parseFloat(ratio) / preRatio;
  let preViewBoxList = this.config.viewBox.split(" ");
  let curWidth = parseInt(parseInt(preViewBoxList[2]) / changeRatio);
  let curHeight = parseInt(parseInt(preViewBoxList[3]) / changeRatio);
  let curX = Math.round(
    parseInt(preViewBoxList[0]) + (parseInt(preViewBoxList[2]) - curWidth) / 2
  );
  let curY = Math.round(
    parseInt(preViewBoxList[1]) + (parseInt(preViewBoxList[3]) - curHeight) / 2
  );
  let curViewBox = `${curX} ${curY} ${curWidth} ${curHeight}`;
  this.nodes.setAttribute("viewBox", curViewBox);
  this.config.viewBox = curViewBox;
  this.config.canvasWidth = curWidth;
  this.config.canvasHeight = curHeight;
  this.noteNode.textContent = `Canvas: ${curWidth} × ${curHeight}`;
  // "Canvas: " + this.config.canvasWidth + " × " + this.config.canvasHeight;
};

/**
 * Bind to zoom in or out events while the mouse is scrolling/wheeling in the SVG area
 */
IbsCharts.prototype.zoomByMouseWheel = function(editable = true) {
  if (editable) {
    if ("object" === typeof this.nodes.onmousewheel) {
      // chrome
      this.nodes.addEventListener("mousewheel", e => {
        e.preventDefault();
        let curRatio = this.config.zoomRatio;
        let ratio = e.wheelDelta > 0 ? curRatio + 0.1 : curRatio - 0.1;
        ratio =
          ratio < this.config.zoomRange[0] ? this.config.zoomRange[0] : ratio;
        ratio =
          ratio > this.config.zoomRange[1] ? this.config.zoomRange[1] : ratio;
        this.zoom(ratio);
      });
    } else {
      // firefox: typeof element.onmousewheel === "undefine"
      this.nodes.addEventListener("DOMMouseScroll", e => {
        e.preventDefault();
        let curRatio = this.config.zoomRatio;
        let ratio = e.detail < 0 ? curRatio + 0.1 : curRatio - 0.1;
        ratio =
          ratio < this.config.zoomRange[0] ? this.config.zoomRange[0] : ratio;
        ratio =
          ratio > this.config.zoomRange[1] ? this.config.zoomRange[1] : ratio;
        this.zoom(ratio);
      });
    }
  }
};

/**
 * Get the gradient url
 * @param {String} type
 * @param {String} color
 */
IbsCharts.prototype.getGradientUrl = function(type, color) {
  color = IbsUtils.getHexColor(color);
  let url = "gradient_" + type + "_" + color;
  if (-1 === this.config.gradientIdList.indexOf(url)) {
    this.config.gradientIdList.push(url);
    let gradient = IbsUtils.createGradient(type, color);
    this.defsNode.appendChild(gradient);
  }
  return url;
};

/**
 * Get the texture url
 * @param {String} type
 * @param {String} color
 */
IbsCharts.prototype.getTextureUrl = function(type, color) {
  color = IbsUtils.getHexColor(color);
  if ("" || "none" === type) {
    return "";
  }
  let url = "texture_" + type + "_" + color;
  if (-1 === this.config.textureIdList.indexOf(url)) {
    this.config.textureIdList.push(url);
    let texture = IbsUtils.createTexture(type, color);
    this.defsNode.appendChild(texture);
  }
  return url;
};

/**
 * show or hide the canvasNote;
 * @param {Boolean} isShow
 */
IbsCharts.prototype.showCanvasNote = function(isShow = true) {
  this.config.showCanvasSize = isShow;
  let display = isShow ? "block" : "none";
  // this.noteNode.setAttribute("display", display);
  this.noteNode.style.display = display;
};

IbsCharts.prototype.resize = function() {
  for (let i = 0; i < this.children.length; i++) {
    this.children[i].update(this.children[i].option);
  }
};

/**
 * create an init config object for the new IBS chart and return
 * @param {HTMLElement} dom
 * @param {Boolean} editable
 */
function createIbsConfig(dom, editable) {
  let config = {
    id: "",
    canvasWidth: 300,
    canvasHeight: 300,
    canvasMinWidth: 300,
    canvasMinHeight: 200,
    viewBox: "0 0 300 200",
    showCanvasSize: false,
    showCanvasBorder: false,
    canvasBorder: null,
    minX: 0,
    maxX: 0,
    minY: 0,
    maxY: 0,
    proteinOrNucleotideStartX: 30,
    proteinOrNucleotideEndX: 270,
    // defaultLength: 1000,
    fixedScaleLength: false,
    scaleLength: 1000,
    lengthRatio: 1000,
    // aspectRatio: "4:3",
    curProteinOrNucleotideY: 0,
    editable: editable,
    gradientIdList: [],
    textureIdList: [],
    // curRatio: 1,
    // wheelZoom: editable,
    zoomRatio: 1,
    zoomRange: [0.2, 5],
    isUndo: false,
    undoStack: [],
    redoStack: [],
    selected: null,
    selectedMolecular: null,
    type: "root",
    heightInterval: 30,
    elementsCount: {
      proteinOrNucleotide: 0,
      line: 0,
      marker: 0,
      bracket: 0,
      text: 0
    },
    isHideAll: true
  };
  let divClientWidth = Math.floor(dom.clientWidth);
  if (config.canvasMinWidth > divClientWidth) {
    dom.style.width = config.canvasMinWidth + "px";
    config.canvasWidth = config.canvasMinWidth;
  } else {
    config.canvasWidth = divClientWidth;
  }
  let divClientHeight = Math.floor(dom.clientHeight);
  if (config.canvasMinHeight > divClientHeight) {
    //dom.style.height = config.canvasHeight + "px";
    config.canvasHeight = config.canvasMinHeight;
  } else {
    config.canvasHeight = divClientHeight;
  }
  let svgViewBox = "0 0 " + config.canvasWidth + " " + config.canvasHeight;
  config.viewBox = svgViewBox;
  config["proteinOrNucleotideStartX"] = parseInt(0.1 * config.canvasWidth);
  config["proteinOrNucleotideEndX"] = parseInt(0.9 * config.canvasWidth);
  config.heightInterval = parseInt(0.05 * config.canvasHeight);
  config.curProteinOrNucleotideY = 8 * config.heightInterval;
  let originalLenRatio =
    (config.proteinOrNucleotideEndX - config.proteinOrNucleotideStartX) /
    config.scaleLength;
  config.lengthRatio =
    originalLenRatio < 0.1
      ? Number(originalLenRatio.toPrecision(2))
      : Number(originalLenRatio.toFixed(2));
  config.id = IbsUtils.getRandomId("IBSCharts");
  dom.innerHTML = "";
  return config;
}

/**
 * Creat a note of the canvas size, like: Canvas: 800 × 600
 * @param {{canvasWidth: Number, canvasHeight: Number, showCanvasSize: Boolean}} config
 */
function createCanvasNote(config) {
  let noteId = IbsUtils.getRandomId("IBSCHarts_canvasNote");
  let noteDiv = document.createElement("div");
  noteDiv.id = noteId;
  let display = config.showCanvasSize ? "block" : "none";
  let style = `display:${display};border:1px solid #333333;width:150px;text-align:center;margin:0;height:12px;position:absolute;line-height:12px;font-size:12px;font-family:Arial;font-weight:bold;`;
  noteDiv.style.cssText = style;
  noteDiv.textContent =
    "Canvas: " + config.canvasWidth + " × " + config.canvasHeight;
  return noteDiv;
}

IbsCharts.prototype.setBasicConfig = function(basicConfig = {}) {
  let keys = [
    "canvasWidth",
    "canvasHeight",
    "viewBox",
    "scaleLength",
    "zoomRatio"
  ];
  if (Object.hasOwnProperty.call(basicConfig, "scaleLength")) {
    this.config.fixedScaleLength = true;
  }
  for (let i = 0; i < keys.length; i++) {
    if (Object.hasOwnProperty.call(basicConfig, keys[i])) {
      this.config[keys[i]] = basicConfig[keys[i]];
    }
  }
  this.config["proteinOrNucleotideStartX"] = parseInt(
    0.1 * this.config.canvasWidth
  );
  this.config["proteinOrNucleotideEndX"] = parseInt(
    0.9 * this.config.canvasWidth
  );
  this.config.heightInterval = parseInt(0.05 * this.config.canvasHeight);
  this.config.curProteinOrNucleotideY = 8 * this.config.heightInterval;
  let originalRatio =
    ((this.config.proteinOrNucleotideEndX -
      this.config.proteinOrNucleotideStartX) *
      this.config.zoomRatio) /
    this.config.scaleLength;
  this.config.lengthRatio =
    originalRatio < 0.1
      ? Number(originalRatio.toPrecision(2))
      : Number(originalRatio.toFixed(2));
  // Number((((this.config.proteinOrNucleotideEndX - this.config.proteinOrNucleotideStartX) * this.config.zoomRatio) / this.config.scaleLength).toFixed(4));
  if (this.noteNode) {
    this.noteNode.textContent =
      "Canvas: " + this.config.canvasWidth + " × " + this.config.canvasHeight;
  }
};

/**
 * download the IBSChrtas in svg, png or jpg format
 * @param {String} type
 */
IbsCharts.prototype.download = function(type = "svg", fileName) {
  // this.showCanvasBorder(false);
  IbsUtils.hideAllBorder();
  type = type.toLowerCase();
  if (!fileName) {
    fileName = "IBSCharts_" + IbsUtils.getNow();
  }

  if ("json" === type) {
    let dataList = [];
    let mode = "protein";
    this.children.forEach(item => {
      if ("protein" === item.config.type || "nucleotide" === item.config.type) {
        mode = item.config.type;
        let pOrNChildren = [];
        item.children.forEach(i => {
          let childOption = i.option;
          pOrNChildren.push({ type: i.config.type, option: childOption });
        });
        item.option.style.align = "custom";
        let pOrNOption = {
          type: item.config.type,
          length: item.config.length,
          option: item.option,
          children: pOrNChildren
        };
        dataList.push(pOrNOption);
      } else {
        // marker, line, text,...
        let option = {
          type: item.config.type,
          option: item.option
        };
        dataList.push(option);
      }
    });
    let jsonObj = {
      // Canvas: [this.config.canvasWidth, this.config.canvasHeight],
      mode: mode,
      config: {
        canvasWidth: this.config.canvasWidth,
        canvasHeight: this.config.canvasHeight,
        viewBox: this.config.viewBox,
        scaleLength: this.config.scaleLength,
        zoomRatio: this.config.zoomRatio
      },
      // viewBox: this.config.viewBox,
      data: dataList
    };
    let jsonStr = JSON.stringify(jsonObj);
    IbsUtils.downloadFile(jsonStr, fileName + ".json");
  } else {
    let preConfig = {
      viewBox: this.config.viewBox,
      zoomRatio: this.config.zoomRatio,
      canvasWidth: this.config.canvasWidth,
      canvasHeight: this.config.canvasHeight,
      proteinOrNucleotideStartX: this.config.proteinOrNucleotideStartX,
      proteinOrNucleotideEndX: this.config.proteinOrNucleotideEndX
    };
    this.adjustCanvasSize();
    if ("jpg" === type) {
      let viewBoxList = this.config.viewBox.split(" ");
      let backgroundAttr = {
        x: parseInt(viewBoxList[0]) - 0.5 * this.config.canvasWidth,
        y: parseInt(viewBoxList[1]) - 0.5 * this.config.canvasHeight,
        width: 2 * this.config.canvasWidth,
        height: 2 * this.config.canvasHeight,
        fill: "#FFFFFF"
      };
      let background = IbsUtils.createSvgElement("rect", backgroundAttr);
      IbsUtils.insertAfter(background, this.defsNode);
      IbsUtils.downloadSvg(this.nodes, fileName, type);
      setTimeout(() => {
        this.nodes.removeChild(background);
      }, 500);
    } else {
      IbsUtils.downloadSvg(this.nodes, fileName, type);
    }

    this.nodes.setAttribute("viewBox", preConfig.viewBox);
    for (const key in preConfig) {
      this.config[key] = preConfig[key];
    }
    this.noteNode.textContent =
      "Canvas: " + this.config.canvasWidth + " × " + this.config.canvasHeight;
  }
  // this.showCanvasBorder(true);
};

/**
 * clear the canvas and return a new instance
 */
IbsCharts.prototype.clear = function() {
  return init(this.parent, true);
};

/**
 * get default option of protein, nucletide or marker
 * @param {String} type
 */
IbsCharts.prototype.getDefaultOption = function(type) {
  // let id = IbsUtils.getRandomId(type);
  let defaultHeight = 25;
  let defaultColor = "#babdb6";
  let defaultEnd = this.config.scaleLength + 1;
  let defaultGradient = "centerToTopAndBottom";
  let defaultOption = {
    id: "",
    displayName: "",
    position: {
      start: {
        site: 1,
        display: "hide"
      },
      end: {
        site: defaultEnd,
        display: "hide"
      }
    },
    coordinate: {
      horizontal: {
        start: this.config.proteinOrNucleotideStartX,
        end: this.config.proteinOrNucleotideEndX,
        isLocked: !this.config.editable
      },
      vertical: {
        start: 0,
        isLocked: !this.config.editable
      }
    },
    style: {
      align: "selfLocation",
      height: defaultHeight,
      fontSize: 12,
      color: defaultColor,
      gradient: defaultGradient,
      texture: {
        type: "none",
        color: "#333333"
      }
    },
    borderStyle: {
      color: "#000000",
      size: 1,
      isDash: false
    }
  };
  switch (type) {
    case "nucleotide":
      defaultOption.id = `Nucleotide${this.config.elementsCount[
        "proteinOrNucleotide"
      ] + 1}`;
      defaultOption.style.height = 5;
      defaultOption.style.color = "#333333";
      defaultOption.style.gradient = "none";
      // defaultOption.position.end.site = 200;
      let temp =
        this.config["curProteinOrNucleotideY"] + this.config.heightInterval;
      if (temp > this.config.canvasHeight) {
        this.config["curProteinOrNucleotideY"] =
          this.config.heightInterval + IbsUtils.getRandomNumber(10, 30);
      } else {
        this.config["curProteinOrNucleotideY"] = temp;
      }
      defaultOption.coordinate.vertical.start = this.config[
        "curProteinOrNucleotideY"
      ];
      // defaultOption.id = IbsUtils.getRandomId("nucleotide");
      break;
    case "marker":
      defaultOption = {
        shape: "circle",
        coordinate: {
          cx: 0.5 * this.config.canvasWidth,
          cy: this.config.curProteinOrNucleotideY + 40
        },
        width: 50,
        height: 50,
        id: "",
        displayName: `Marker${this.config.elementsCount["marker"] + 1}`,
        text: {
          content: "",
          position: "center", //number or string [center,top,bottom]
          style: {
            fontSize: 14,
            fontFamily: "Arial",
            color: "#333333",
            borderColor: "black",
            fontStyle: "normal", //斜体
            fontWeight: "normal", //粗体
            // textDecoration: "none",
            // horizontalAlign: "middle", //textAnchor
            // verticalAlign: "center", //
            rotate: 0
          }
        },
        style: {
          color: "#bc0000",
          borderColor: "#000000",
          borderSize: 1,
          isDash: false,
          rotate: 0,
          gradient: "centerToAllAround",
          texture: {
            type: "none",
            color: "#4389b4"
          }
        },
        editable: this.config.editable
      };
      break;
    case "line":
      let defaultY = this.config["curProteinOrNucleotideY"];
      // let defaultLineLength = 0.15 * this.config.canvasWidth;
      defaultOption = {
        type: "straight", //straight, polyline, curve
        // position: "200,200 400,0",
        position: `${0.4 * this.config.canvasWidth},${defaultY} ${0.6 *
          this.config.canvasWidth},${defaultY}`,
        pointsNum: 2,
        id: "",
        displayName: `Line${this.config.elementsCount["line"] + 1}`,
        style: {
          isDash: false,
          left: "none",
          right: "none",
          color: "#000000",
          size: 1
        },
        editable: this.config.editable
      };
      break;
    case "text":
      defaultOption = {
        content: "Text",
        position: {
          x: 0.5 * this.config.canvasWidth,
          y: this.config["curProteinOrNucleotideY"] + 30
        },
        displayName: `Text${this.config.elementsCount["text"] + 1}`,
        id: "",
        style: {
          fontSize: 18,
          fontFamily: "Arial",
          color: "#000000",
          fontStyle: "normal",
          fontWeight: "normal",
          textDecoration: "none", //textDecoration
          horizontalAlign: "middle", // textAnchor
          verticalAlign: "center", //[center, buttom, top]
          rotate: 0
        },
        borderStyle: {
          //外部矩形边框
          color: "#333333", //stroke
          size: 0, //strokeWidth
          isDash: false //dashArray=1,1
        },
        editable: this.config.editable
      };
      break;
    case "bracket":
      defaultOption = {
        id: "",
        displayName: `Bracket${this.config.elementsCount["bracket"] + 1}`,
        shape: "brace",
        coordinate: {
          cx: 0.5 * this.config.canvasWidth,
          cy: this.config["curProteinOrNucleotideY"] + 30
        },
        style: {
          color: "#000000",
          size: 1,
          rotate: 0,
          width: 100
        },
        editable: this.config.editable
      };
      break;
    default:
      // protein
      defaultOption.id = `Protein${this.config.elementsCount[
        "proteinOrNucleotide"
      ] + 1}`;
      let tempHeight =
        this.config["curProteinOrNucleotideY"] + this.config.heightInterval;
      if (tempHeight > this.config.canvasHeight) {
        this.config["curProteinOrNucleotideY"] =
          this.config.heightInterval + IbsUtils.getRandomNumber(10, 30);
      } else {
        this.config["curProteinOrNucleotideY"] = tempHeight;
      }
      // this.config["curProteinOrNucleotideY"] += parseInt(this.config.canvasHeight / 10);
      defaultOption.coordinate.vertical.start = this.config[
        "curProteinOrNucleotideY"
      ];
      // defaultOption.id = IbsUtils.getRandomId("protein");
      break;
  }
  return defaultOption;
};

/**
 * Create a basic protein or nucleotide chart, and you can add other components to it
 * @param {{id: String, position: {start: {startSite: *, display: String}, end: {endSite: *, display: String} },
 * coordinate: {horizontal: {start: String, end: String, isLocked: Boolean,}, vertical: {start: String, isLocked: Boolean}, },
 * style: {align: String, height: Number, fontSize: Number, color: String, gradient:String, texture: {type: String, color: String}},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} proteinOption
 */
IbsCharts.prototype.createProteinOrNucleotide = function(
  option,
  type = "protein"
) {
  type =
    !type || ("protein" !== type && "nucleotide" !== type) ? "protein" : type;
  let defaultOption = this.getDefaultOption(type);
  option = IbsUtils.supplementToStandardObj(defaultOption, option);
  option.position.start.site = parseInt(option.position.start.site);
  option.position.end.site = parseInt(option.position.end.site);
  option.position.start.site = Math.min(
    option.position.start.site,
    option.position.end.site
  );
  option.position.end.site = Math.max(
    option.position.start.site,
    option.position.end.site
  );
  this.config.elementsCount.proteinOrNucleotide += 1;
  let length = option.position.end.site - option.position.start.site;
  if (
    !this.config.fixedScaleLength &&
    1 === this.config.elementsCount.proteinOrNucleotide &&
    length !== this.config.scaleLength
  ) {
    this.config.scaleLength = length;
    let originalLenRatio =
      (this.config.proteinOrNucleotideEndX -
        this.config.proteinOrNucleotideStartX) /
      this.config.scaleLength;
    this.config.lengthRatio =
      originalLenRatio < 0.1
        ? Number(originalLenRatio.toPrecision(2))
        : Number(originalLenRatio.toFixed(2));
    // Number(((this.config.proteinOrNucleotideEndX - this.config.proteinOrNucleotideStartX) / this.config.scaleLength).toFixed(4));
  }
  let pOrN;
  switch (type) {
    case "nucleotide":
      pOrN = new NucleotideModel.Nucleotide(option, this);
      break;
    case "protein":
      pOrN = new ProteinModel.Protein(option, this);
      break;
    default:
      break;
  }
  this.nodes.appendChild(pOrN.nodes);
  this.children.push(pOrN);
  let horizontalDraggable =
    false === option.coordinate.horizontal.isLocked ? true : false;
  let verticalDraggable =
    false === option.coordinate.vertical.isLocked ? true : false;
  pOrN.draggable({
    horizontal: horizontalDraggable,
    vertical: verticalDraggable
  });
  this.config.isUndo = false;
  let operation = {
    target: pOrN,
    cmd: "create",
    args: JSON.stringify(pOrN.option),
    nodesIndex: pOrN.nodesIndex
  };
  this.config.undoStack.push(operation);
  return pOrN;
};

/**
 * create a text on the canvas
 * @param {*} option
 */
IbsCharts.prototype.createText = function(option) {
  option.editable = this.config.editable;
  let defaultOption = this.getDefaultOption("text");
  let textOption = IbsUtils.supplementToStandardObj(defaultOption, option);
  let text = new Text(textOption, this);
  this.children.push(text);
  this.config.elementsCount["text"] += 1;
  this.nodes.appendChild(text.nodes);
  this.config.isUndo = false;
  let operation = {
    target: text,
    cmd: "create",
    args: JSON.stringify(text.option),
    nodesIndex: text.nodesIndex
  };
  this.config.undoStack.push(operation);
  return text;
};

IbsCharts.prototype.createMarker = function(option) {
  let marker = new Marker(option, this);
  this.children.push(marker);
  this.config.elementsCount.marker += 1;
  this.nodes.appendChild(marker.nodes);
  this.config.isUndo = false;
  let operation = {
    target: marker,
    cmd: "create",
    args: JSON.stringify(marker.option),
    nodesIndex: marker.nodesIndex
  };
  this.config.undoStack.push(operation);
  return marker;
};

IbsCharts.prototype.createLine = function(option) {
  // option.editable = true;
  let defaultOption = this.getDefaultOption("line");
  option = IbsUtils.supplementToStandardObj(defaultOption, option);
  let line = new Line(option, this);
  this.children.push(line);
  this.config.elementsCount.line += 1;
  this.nodes.appendChild(line.nodes);
  this.config.isUndo = false;
  let operation = {
    target: line,
    cmd: "create",
    args: JSON.stringify(line.option),
    nodesIndex: line.nodesIndex
  };
  this.config.undoStack.push(operation);
  return line;
};

IbsCharts.prototype.createBracket = function(option) {
  // option.editable = true;
  let defaultoption = this.getDefaultOption("bracket");
  let bracketOption = IbsUtils.supplementToStandardObj(defaultoption, option);
  let bracket = new Bracket(bracketOption, this);
  this.children.push(bracket);
  this.config.elementsCount.bracket += 1;
  this.nodes.appendChild(bracket.nodes);
  this.config.isUndo = false;
  let operation = {
    target: bracket,
    cmd: "create",
    args: JSON.stringify(bracket.option),
    nodesIndex: bracket.nodesIndex
  };
  this.config.undoStack.push(operation);
  return bracket;
};

/**
 * Create a basic protein chart, and you can add other components to it
 * @param {{id: String, displayName: String, position: {start: {startSite: *, display: String}, end: {endSite: *, display: String} },
 * coordinate: {horizontal: {start: String, end: String, isLocked: Boolean,}, vertical: {start: String, isLocked: Boolean}, },
 * style: {align: String, height: Number, fontSize: Number, color: String, gradient:String, texture: {type: String, color: String}},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 */
IbsCharts.prototype.createProtein = function(option) {
  return this.createProteinOrNucleotide(option, "protein");
};

/**
 * Create a basic nucleotide chart, and you can add other components to it
 * @param {{id: String, displayName: String, position: {start: {startSite: *, display: String}, end: {endSite: *, display: String} },
 * coordinate: {horizontal: {start: String, end: String, isLocked: Boolean,}, vertical: {start: String, isLocked: Boolean}, },
 * style: {align: String, height: Number, fontSize: Number, color: String, gradient:String, texture: {type: String, color: String}},
 * borderStyle: {color: String, size: Number, isDash: Boolean}}} option
 */
IbsCharts.prototype.createNucleotide = function(option) {
  return this.createProteinOrNucleotide(option, "nucleotide");
};

/**
 * Undo the previous step
 */
IbsCharts.prototype.undo = function() {
  let operation = this.config.undoStack.pop();
  if (operation) {
    switch (operation.cmd) {
      case "delete":
        let parent = operation.target.parent;
        if (parent.children[operation.nodesIndex]) {
          parent.nodes.insertBefore(
            operation.target.nodes,
            parent.children[operation.nodesIndex].nodes
          );
        } else {
          parent.nodes.appendChild(operation.target.nodes);
        }
        parent.children.splice(operation.nodesIndex, 0, operation.target);
        // this.config.undoStack.pop();
        break;
      case "create":
        operation.target.delete();
        this.config.undoStack.pop();
        break;
      default:
        // other interaction: update, drag, ...
        let option = JSON.parse(operation.args);
        operation.target.update(option);
        operation = this.config.undoStack.pop();
        break;
    }
    this.config.redoStack.push(operation);
    this.config.isUndo = true;
  } else {
    console.log("can't undo more!");
  }
  // console.log("undoStack: ");
  // console.log(this.config.undoStack);
  // console.log("redoStack: ");
  // console.log(this.config.redoStack);
};

/**
 * Redo the operation that was undone in the previous step
 */
IbsCharts.prototype.redo = function() {
  if (this.config.isUndo) {
    let operation = this.config.redoStack.pop();
    if (operation) {
      switch (operation.cmd) {
        case "delete":
          operation.target.delete();
          break;
        case "create":
          let parent = operation.target.parent;
          if (parent.children[operation.nodesIndex]) {
            parent.nodes.insertBefore(
              operation.target.nodes,
              parent.children[operation.nodesIndex].nodes
            );
          } else {
            parent.nodes.appendChild(operation.target.nodes);
          }
          parent.children.splice(operation.nodesIndex, 0, operation.target);
          this.config.undoStack.push(operation);
          break;
        default:
          // other interaction: update, drag, ...
          operation.target.update(JSON.parse(operation.args));
          break;
      }
    } else {
      console.log("can't redo more!");
    }
    this.config.isUndo = true;
  } else {
    this.config.redoStack.length = 0; // clear the redoStack
    console.log("redo stack is empty!");
  }
};

/**
 * Hide the element's selected border by className: nodeBorder
 */
IbsCharts.prototype.hideSelectedBorder = function() {
  document.addEventListener("mousedown", e => {
    IbsUtils.hideAllBorder();
    this.config.isHideAll = true;
  });
};

/**
 * Create a Ibs element
 * @param {String} type
 * @param {Object} option
 */
IbsCharts.prototype.createIbsElement = function(type, option) {
  switch (type) {
    case "protein":
      this.createProtein(option);
      break;
    case "nucleotide":
      this.createNucleotide(option);
      break;
    case "text":
      this.createText(option);
      break;
    case "marker":
      this.createMarker(option);
      break;
    case "line":
      this.createLine(option);
      break;
    case "bracket":
      this.createBracket(option);
      break;
    case "domain":
      this.config.selectedMolecular.createDomain(option);
      break;
    case "cutline":
      this.config.selectedMolecular.createCutline(option);
      break;
    case "site":
      this.config.selectedMolecular.createSite(option);
      break;
    default:
      break;
  }
};

/**
 * crea
 * @param {Object} jsonObj
 */
IbsCharts.prototype.loadJson = function(jsonObj) {
  if (Object.hasOwnProperty.call(jsonObj, "config")) {
    this.setBasicConfig(jsonObj.config);
  }
  jsonObj.data.forEach(item => {
    let defaultOption = this.getDefaultOption(item.type);
    item.option = IbsUtils.supplementToStandardObj(defaultOption, item.option);
    this.createIbsElement(item.type, item.option);
    if (item.children && 0 < item.children.length) {
      item.children.forEach(subItem => {
        try {
          if ("site" == subItem.type) {
            let subItemDefaultOption = this.getDefaultOption(subItem.type);
            subItem.option = IbsUtils.supplementToStandardObj(
              subItemDefaultOption,
              subItem.option
            );
          }
          this.createIbsElement(subItem.type, subItem.option);
        } catch (error) {
          console.log(error);
        }
      });
    }
  });
  this.centerCanvas();
  // this.adjustCanvasSize();
  IbsUtils.hideAllBorder();
};
export default { init };