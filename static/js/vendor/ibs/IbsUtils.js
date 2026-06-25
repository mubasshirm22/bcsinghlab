/*
 * Common methods component of IBS
 */
import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";

/**
 * Automatically adds default values for properties that were not set in the original object
 * @param {{}} standrdObj 
 * @param {{}} originalObj 
 */
function supplementToStandardObj(standrdObj, originalObj) {
  if (standrdObj && typeof standrdObj === "object") {
    if (originalObj && typeof originalObj === "object") {
      for (const key in standrdObj) {
        if (standrdObj.hasOwnProperty(key)) {
          if (originalObj.hasOwnProperty(key)) {
            // Determine whether the obj child element is an object, if it is, recursively detection and set default value
            if (standrdObj[key] && typeof standrdObj[key] === "object") {
              supplementToStandardObj(standrdObj[key], originalObj[key])
            }
          } else {
            // If not, simply add the property and the default value
            // originalObj[key] = standrdObj[key];
            originalObj[key] = deepCopyObj(standrdObj[key]);
          }
        }
      }
      return originalObj;
    } else {
      return deepCopyObj(standrdObj);
    }
  } else {
    return false;
  }
}

/** 
 * Deep copy an object or other simple value
 * @param {*} obj 
 */
function deepCopyObj(obj) {
  let objClone;
  if (obj && typeof obj === "object") {
    objClone = Array.isArray(obj) ? [] : {};
    for (const key in obj) {
      if (obj.hasOwnProperty(key)) {
        objClone[key] = deepCopyObj(obj[key]);
      }
    }
  } else {
    objClone = obj;
  }
  return objClone;
}

/**
 * get the first child element node
 * @param {HTMLElement} element 
 */
function getFirstElement(element) {
  // If the browser supports using "firstElementChild" to get the first tag node
  if (element.firstElementChild) {
    return element.firstElementChild; // Returns the retrieved tag node object
  } else { // If not, use "firstChild" to get the firstChild of the current tag
    var ele = element.firstChild; // To get the first child node
    while (ele && ele.nodeType !== 1) { // The first ele is to make sure that there is this object, and the second is that the current node type is not a element node
      ele = ele.nextSibling; // If not, proceed to get the next node
    }
    // Returns the first child tag node retrieved
    return ele;
  }
}

/**
 * To prevent the bubbling up
 * @param {MouseEvent} e 
 */
function stopBubble(e) {
  // To stop bubbling: 1. For non-IE: stopPropagation(). 2. For Internet Explorer: Set cancelBubble property to true
  e.stopPropagation ? e.stopPropagation() : e.cancelBubble = true;
}

/**
 * Get a random number between min and Max
 * @param {Number} min 
 * @param {Number} max 
 */
function getRandomNumber(min = 0, max = 1) {
  let n = Math.random() * (max - min);
  return Math.floor(n + min);
}

/**
 * Generates a random ID with the specified prefix: "finger"
 * @param {String} finger 
 */
function getRandomId(finger = "") {
  let id = finger + "_";
  const CHARTS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";
  for (let i = 0; i < 5; i++) {
    let n = getRandomNumber(0, CHARTS.length - 1);
    id += CHARTS[n];
  }
  // let d = new Date();
  // let milliseconds = d.getTime().toString().slice(-3);
  // id += "_" + milliseconds;
  return id;
}

/**
 * create a SVG element and set atrributes
 * @param {String} tag 
 * @param {{}} attrObj 
 * @return {SVGAElement}
 */
function createSvgElement(tag, attrObj = {}) {
  var svgTag = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const key in attrObj) {
    // Cycle setting attributes
    svgTag.setAttribute(key, attrObj[key]);
  }
  return svgTag;
}

/**
 * show the label with the display: "top" or "bottom"
 * @param {{id: String, type:String, position:String, display: String, label: String, x: Number, y: Number, lineLength: Number}} labelOption 
 */
function createStartOrEndLabel(labelOption) {
  let visibility = "hide" === labelOption.display ? "hidden" : "visible";
  let labelGAttr = {
    id: labelOption.id,
    type: "label",
    visibility: visibility
  };
  let labelGTag = createSvgElement("g", labelGAttr);
  let labelAttrs = getLabelAttr(labelOption);
  let lineAttr = {
    d: labelAttrs.d,
    style: `stroke:#000000;fill:none;`
  };
  let line = createSvgElement("path", lineAttr);
  labelGTag.appendChild(line);
  let labelText = createSvgElement("text", labelAttrs.labelTextAttr);
  labelText.textContent = labelOption.label;
  labelGTag.appendChild(labelText);
  return labelGTag;
}

/**
 * @param {{id: String, type:String, position:String, display: String, label: String, x: Number, y: Number, lineLength: Number}} labelOption 
 */
function getLabelAttr(labelOption) {
  let vDirection = "bottom" === labelOption.display ? "" : "-";
  let d = `M${labelOption.x},${labelOption.y} v${vDirection}${labelOption.lineLength}`;
  let lineEndY = "bottom" === labelOption.display ? labelOption.y + labelOption.lineLength : labelOption.y - labelOption.lineLength;
  let dy = "bottom" === labelOption.display ? "1em" : "-0.2em";
  let labelTextX = labelOption.x;
  if ("polyline" === labelOption.type) {
    let xDirection = "-";
    let yDirection = "";
    if ("end" === labelOption.position) {
      xDirection = "";
      labelTextX += labelOption.lineLength;
    } else {
      labelTextX -= labelOption.lineLength;
    }
    if ("bottom" === labelOption.display) {
      lineEndY += labelOption.lineLength;
    } else {
      lineEndY -= labelOption.lineLength;
      yDirection = "-"
    }
    d += ` l${xDirection}${labelOption.lineLength},${yDirection}${labelOption.lineLength}`;
  }
  let labelTextAttr = {
    id: labelOption.id + "_label",
    x: labelTextX,
    y: lineEndY,
    dy: dy,
    "text-anchor": "middle"
  };
  return {
    d: d,
    labelTextAttr: labelTextAttr
  }
}

/**
 * 
 * @param {{id: String, type:String, position:String, display: String, label: String, x: Number, y: Number, lineLength: Number}} labelOption 
 * @param {SVGAElement} parent 
 */
function updateLabel(labelOption, parent) {
  let label = document.getElementById(labelOption.id);
  if (label) { // if exists
    if ("hide" === labelOption.display) {
      label.setAttribute("visibility", "hidden");
    } else {
      label.setAttribute("visibility", "visible");
      let children = label.children;
      let labelAttrs = getLabelAttr(labelOption);
      children[0].setAttribute("d", labelAttrs.d);
      children[1].setAttribute("x", labelAttrs.labelTextAttr.x);
      children[1].setAttribute("y", labelAttrs.labelTextAttr.y);
      children[1].setAttribute("dy", labelAttrs.labelTextAttr.dy);
      children[1].textContent = labelOption.label;
    }
  } else { // if not exists
    label = createStartOrEndLabel(labelOption);
    parent.appendChild(label);
  }
  return label;
}

/**
 * 
 * @param {*} domain 
 * @param {ProteinOrNucleotide} parent 
 */
function getDomainLabelOption(domain, parent) {
  let labelLength = Math.max(6, parent.option.style.height / 2);
  let PY = parent.option.coordinate.vertical.start;
  let PHeight = parent.option.style.height;
  let startLabelOption = {
    id: domain.config.groupId + "_start",
    type: "line",
    position: "start",
    display: domain.option.position.start.display,
    label: domain.option.position.start.site,
    x: domain.config.startX,
    y: "bottom" === domain.option.position.start.display ? PY + PHeight : PY,
    lineLength: labelLength
  };
  let endLabelOption = {
    id: domain.config.groupId + "_end",
    type: "line",
    position: "end",
    display: domain.option.position.end.display,
    label: domain.option.position.end.site,
    x: domain.config.startX + domain.config.width,
    y: "bottom" === domain.option.position.end.display ? PY + PHeight : PY,
    lineLength: labelLength
  };
  return {
    startLabelOption: startLabelOption,
    endLabelOption: endLabelOption
  };
}

/**
 * Generates a boundary consisting of four vertices
 * @param {[Number, Number][]} coordinates 
 */
function createBorder(coordinates) {
  let borderGroupAttr = {
    class: "nodeBorder",
    visibility: "hidden"
  };
  let borderGroup = createSvgElement("g", borderGroupAttr);
  coordinates.forEach(coordinate => {
    let circleAttr = {
      r: 4,
      cx: coordinate[0],
      cy: coordinate[1],
      fill: "#ffa500"
    };
    let circle = createSvgElement("circle", circleAttr);
    borderGroup.appendChild(circle);
  });
  return borderGroup;
}

/**
 * Show the border of the target element by mousedown, and hide others
 * @param {SVGAElement} targetNode 
 * @param {SVGAElement} borderNode 
 */
function showBorderByMousedown(targetNode, borderNode) {
  targetNode.addEventListener("mousedown", (e) => {
    stopBubble(e);
    showBorderAlone(borderNode);
  })
}

/**
 * Show the border of the target element and hide others
 * @param {SVGAElement} borderNode 
 */
function showBorderAlone(borderNode) {
  borderNode.setAttribute("visibility", "visible");
  let nodeBorders = document.getElementsByClassName("nodeBorder");
  for (let i = 0; i < nodeBorders.length; i++) {
    if (borderNode !== nodeBorders[i]) {
      nodeBorders[i].setAttribute("visibility", "hidden");
    }
  }
}

/**
 * 
 * @param {{parent:{},nodes:SVGAElement, config: {}, option: {}}} target 
 */
function deleteObj(target) {
  // console.log("delete: ");
  let index = target.parent.children.indexOf(target);
  if (-1 !== index) {
    target.parent.children.splice(index, 1);
    // target.nodes = null;
    target.parent.nodes.removeChild(target.nodes);
    // for (const key in target) {
    //   if (target.hasOwnProperty(key)) {
    //     target[key] = null;
    //   }
    // }
  } else {
    console.log("this object is null!");
  }
  return target;
}

/**
 * Get the hexadecimal color code from rgb code
 * @param {Number} r 
 * @param {Number} g 
 * @param {Number} b 
 */
function rgbToHex(r, g, b) {
  r = Math.min(255, Math.max(0, r));
  g = Math.min(255, Math.max(0, g));
  b = Math.min(255, Math.max(0, b));
  let result = ((r << 16) | (g << 8) | b).toString(16);
  let prefix = "#";
  let dif = 6 - result.length;
  for (let i = 0; i < dif; i++) {
    prefix += "0";
  }
  return (prefix + result);
}

/**
 * Get the hexadecimal color code
 * @param {String} color 
 */
function getHexColor(color) {
  color = color.replace(/\s/g, "").toLowerCase();
  if ("#" === color.charAt(0) || "rgb" !== color.substr(0, 3)) {
    return color;
  } else {
    let rgb = color.split(/\D+/);
    let rgbList = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
      let index = i + 1;
      if (rgb[index]) {
        rgbList[i] = parseInt(rgb[index]);
        // rgbList[i] = Math.max(0, rgbList[i]);
        // rgbList[i] = Math.min(255, rgbList[i]);
      }
    }
    return rgbToHex(...rgbList);
  }
}

function createGradient(type, color) {
  let tag = ("centerToAllAround" === type || "allAroundToCenter" === type) ? "radialGradient" : "linearGradient";
  let gradientAttr = {
    id: "gradient_" + type + "_" + color,
    x1: "0%",
    y1: "0%",
    x2: "0%",
    y2: "100%",
    spreadMethod: "pad"
  };
  let childrenAttrList = [{
    "offset": "0%",
    "stop-color": "#ffffff"
  }, {
    "offset": "100%",
    "stop-color": color
  }];
  switch (type) {
    case "topAndBottomToCenter":
      childrenAttrList[1] = {
        "offset": "50%",
        "stop-color": color
      };
      childrenAttrList[2] = {
        "offset": "100%",
        "stop-color": "#ffffff"
      };
      break;
    case "centerToLeftAndRight":
      gradientAttr.x2 = "100%";
      gradientAttr.y2 = "0%";
      childrenAttrList[0] = {
        "offset": "0%",
        "stop-color": color
      };
      childrenAttrList[1] = {
        "offset": "50%",
        "stop-color": "#ffffff"
      };
      childrenAttrList[2] = {
        "offset": "100%",
        "stop-color": color
      };
      break;
    case "leftAndRightToCenter":
      gradientAttr.x2 = "100%";
      gradientAttr.y2 = "0%";
      childrenAttrList[1] = {
        "offset": "50%",
        "stop-color": color
      };
      childrenAttrList[2] = {
        "offset": "100%",
        "stop-color": "#ffffff"
      };
      break;
    case "centerToAllAround":
      gradientAttr = {
        id: "gradient_" + type + "_" + color,
        r: "50%"
      };
      break;
    case "allAroundToCenter":
      gradientAttr = {
        id: "gradient_" + type + "_" + color,
        r: "50%"
      };
      childrenAttrList[0] = {
        "offset": "0%",
        "stop-color": color
      };
      childrenAttrList[1] = {
        "offset": "100%",
        "stop-color": "#ffffff"
      };
      break;
    case "topToBottom":
      break;
    case "bottomToTop":
      childrenAttrList[0] = {
        "offset": "0%",
        "stop-color": color
      };
      childrenAttrList[1] = {
        "offset": "100%",
        "stop-color": "#ffffff"
      };
      break;
    case "topLeftToBottomRight":
      gradientAttr.x2 = "100%";
      break;
    case "topRightToBottomLeft":
      gradientAttr.x1 = "100%";
      break;
    case "bottomLeftToTopRight":
      gradientAttr.y1 = "100%";
      gradientAttr.x2 = "100%";
      gradientAttr.y2 = "0%";
      break;
    case "bottomRightToTopLeft":
      gradientAttr.x1 = "100%";
      gradientAttr.y1 = "100%";
      gradientAttr.y2 = "0%";
      break;
    case "leftToRight":
      gradientAttr.x2 = "100%";
      gradientAttr.y2 = "0%";
      break;
    case "rightToLeft":
      gradientAttr.x1 = "100%";
      gradientAttr.y2 = "0%";
      break;
    case "centerToTopAndBottom":
      childrenAttrList[0] = {
        "offset": "0%",
        "stop-color": color
      };
      childrenAttrList[1] = {
        "offset": "50%",
        "stop-color": "#ffffff"
      };
      childrenAttrList[2] = {
        "offset": "100%",
        "stop-color": color
      };
      break;
    default: //  "" || "none"
      childrenAttrList.splice(0, 1);
      break;
  }
  let gradient = createSvgElement(tag, gradientAttr);
  childrenAttrList.forEach(e => {
    gradient.appendChild(createSvgElement("stop", e));
  });
  return gradient;
}

function createTexture(type, color) {
  if ("crystal" === type) {
    return createCrystalTexture(color);
  }
  let patternAttr = {
    id: "texture_" + type + "_" + color,
    x: "0",
    y: "0",
    width: "5",
    height: "5",
    patternUnits: "userSpaceOnUse"
  };
  let pathAttr = {
    d: "M0,0 h5",
    stroke: color
  };
  switch (type) {
    case "horizontalLine":
      break;
    case "verticalLine":
      pathAttr.d = "M0,0 v5";
      break;
    case "horizontalDashLine":
      pathAttr.d = "M0,0 h2.5";
      break;
    case "verticalDashLine":
      pathAttr.d = "M0,0 v2.5";
      break;
    case "leftDiagonal":
      patternAttr["patternTransform"] = "rotate(135)";
      break;
    case "rightDiagonal":
      patternAttr["patternTransform"] = "rotate(45)";
      break;
    case "leftDashDiagonal":
      pathAttr.d = "M0,0 h2.5";
      patternAttr["patternTransform"] = "rotate(135)";
      break;
    case "rightDashDiagonal":
      pathAttr.d = "M0,0 h2.5";
      patternAttr["patternTransform"] = "rotate(45)";
      break;
    default: // "" || "none"
      pathAttr.d = "";
      break;
  }
  let path = createSvgElement("path", pathAttr);
  let pattern = createSvgElement("pattern", patternAttr);
  pattern.appendChild(path);
  return pattern;
}

/**
 * 
 * @param {String} color 
 */
function createCrystalTexture(color) {
  let gradientAttr = {
    id: "texture_crystal_" + color,
    x1: "0%",
    y1: "0%",
    x2: "0%",
    y2: "100%",
    spreadMethod: "pad"
  };
  let r = parseInt(color.slice(1, 3), 16);
  let g = parseInt(color.slice(3, 5), 16);
  let b = parseInt(color.slice(5), 16);
  let color1 = rgbToHex(r + 155, g + 155, b + 155);
  let color2 = rgbToHex(r + 50, g + 50, b + 50);
  let color3 = rgbToHex(r + 38, g + 38, b + 38);
  let childrenAttrList = [{
      "offset": "0%",
      "stop-color": color1
    },
    {
      "offset": "49%",
      "stop-color": color2
    },
    {
      "offset": "50%",
      "stop-color": color
    },
    {
      "offset": "100%",
      "stop-color": color3
    }
  ];
  let gradient = createSvgElement("linearGradient", gradientAttr);
  childrenAttrList.forEach(i => {
    gradient.appendChild(createSvgElement("stop", i));
  });
  return gradient;
}

function layerMove(target, direction = "up") {
  let index = target.parent.children.indexOf(target);
  if (-1 !== index) {
    if ("down" === direction) {
      if (0 < index) {
        let lastObj = target.parent.children[index - 1];
        target.parent.nodes.insertBefore(target.nodes, lastObj.nodes);
        target.parent.children.splice(index, 1);
        target.parent.children.splice(index - 1, 0, target);
      }
    } else if ("up" === direction) { 
      if (target.parent.children.length - 1 > index) {
        let nextObj = target.parent.children[index + 1];
        insertAfter(target.nodes, nextObj.nodes);
        target.parent.children.splice(index, 1);
        target.parent.children.splice(index + 1, 0, target);
      }
    } else if ("back" === direction) {
      if (0 < index) {
        target.parent.nodes.insertBefore(target.nodes, target.parent.children[0].nodes);
        target.parent.children.splice(index, 1);
        target.parent.children.splice(0, 0, target);
      }
    } else if ("front" === direction) {
      if (target.parent.children.length - 1 > index) {
        target.parent.nodes.removeChild(target.nodes);
        target.parent.nodes.appendChild(target.nodes);
        target.parent.children.splice(index, 1);
        target.parent.children.push(target);
      }
    }
  }
}

/**
 * calculate the angel formed by 3 points
 * @param {Number} x1 the index1 control point x coordinate
 * @param {Number} y1 the index1 control point y coordinate
 * @param {Number} x2 the mouseEvent box.x
 * @param {Number} y2 the mouseEvent box.y
 * @param {Number} x0 the angel point x coordinate
 * @param {Number} y0 the angel point y coordinate
 */
function Cosines(x1, y1, x2, y2, x0, y0) { //有时会有bug，还在找
  var distance1 = Math.sqrt(Math.pow(x1 - x0, 2) + Math.pow(y1 - y0, 2));
  var distance2 = Math.sqrt(Math.pow(x2 - x0, 2) + Math.pow(y2 - y0, 2));
  var dotProduct = (x1 - x0) * (x2 - x0) + (y1 - y0) * (y2 - y0);
  var Rad = Math.acos(dotProduct / distance1 / distance2);
  var flag = x2 >= x0 ? 0 : 1; //还要考虑
  var Angel = Rad / (2 * Math.PI) * 360;
  if (flag) {
    Angel = 0 - Angel;
  }
  // var flag = (typeof value) === 'number' && window.isNaN(value);
  // if (flag) {
  //   console.log(x1, x2, x0, y1, y2, y0);
  //   console.log("error!!!!!");
  // }
  return Angel;
}
/**
 * get the coordinate box of the target node in the reference node coordinate system
 * @param {HTMLElement} target 
 * @param {HTMLElement} reference 
 */
function calAbsoluteCoordinate(target, reference) {
  var box = SVG(target).rbox(SVG(reference));
  return box;
}
/**
 * add mousedown event that display the target node 
 * @param {HTMLElement} targetNode 
 * @param {HTMLElement} clickNode 
 */
function HideAndDisplayControlBox(targetNode, clickNode, parentNode) {
  //var svg = document.querySelector("svg");
  //console.log(svg);
  targetNode.style.visibility = "hidden";
  clickNode.addEventListener("mousedown", function (e) {
    targetNode.style.visibility = "visible";
  }, true)
  parentNode.addEventListener("mousedown", function (e) {
    targetNode.style.visibility = "hidden";
  }, true)
}
/**
 * calculate 8 points coordinate on the sides and angel of a rectangle
 * @param {Number} cx center x coordinate
 * @param {Number} cy center y coordinate
 * @param {Number} width 
 * @param {Number} height 
 */
function calRectEightPoints(cx, cy, width, height) {
  let order = [0, 1, 2, 5, 8, 7, 6, 3];
  var sidePoints = [];
  var centerPoint = [];
  let initX = cx - width * 0.5;
  let initY = cy - height * 0.5;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let num = 3 * i + j;
      let x = initX + width / 2 * j;
      let y = initY + height / 2 * i;
      if (num == 4) {
        centerPoint[0] = x;
        centerPoint[1] = y
        continue;
      }
      let index = order.indexOf(num);
      sidePoints[index] = [];
      sidePoints[index][0] = x;
      sidePoints[index][1] = y;
    }
  }
  return {
    sidePoints,
    centerPoint
  };
}

function calRectEightPointsWithRotate(cx, cy, width, height, rotate) {
  cx = parseFloat(cx);
  cy = parseFloat(cy);
  let order = [0, 1, 2, 5, 8, 7, 6, 3];
  var sidePoints = [];
  var centerPoint = [];
  let initX = cx - width * 0.5;
  let initY = cy - height * 0.5;
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let num = 3 * i + j;
      let x = initX + width / 2 * j;
      let y = initY + height / 2 * i;
      // let r = Math.sqrt(Math.pow((x-cx),2) + Math.pow((y-cy),2))
      if (num == 4) {
        centerPoint[0] = x;
        centerPoint[1] = y
        continue;
      }
      let index = order.indexOf(num);
      sidePoints[index] = [];
      sidePoints[index][0] = cx + (x - cx) * Math.cos(rotate / 180 * Math.PI) - (y - cy) * Math.sin(rotate / 180 * Math.PI);
      sidePoints[index][1] = cy + (y - cy) * Math.cos(rotate / 180 * Math.PI) + (x - cx) * Math.sin(rotate / 180 * Math.PI);
    }
  }
  return {
    sidePoints,
    centerPoint
  };
}

function calculatePointAfterRotate(cx, cy, x, y, rotate) {
  cx = parseFloat(cx);
  cy = parseFloat(cy);
  x = parseFloat(x);
  y = parseFloat(y);
  let resX = cx + (x - cx) * Math.cos(rotate / 180 * Math.PI) - (y - cy) * Math.sin(rotate / 180 * Math.PI);
  let resY = cy + (y - cy) * Math.cos(rotate / 180 * Math.PI) + (x - cx) * Math.sin(rotate / 180 * Math.PI);
  return [resX, resY];
}

/**
 * transform the ibs style option to svg sytle option
 * @param {{}} styleOption 
 */
function standarizeStyleOption(styleOption) {
  let option = {};
  for (let key in styleOption) {
    let value = styleOption[key];
    switch (key) {
      case "color":
        if (styleOption.hasOwnProperty("right")) {
          option.stroke = value;
        } else {
          option.fill = value;
        }
        break;
      case "borderColor":
        option.stroke = value;
        break;
      case "borderWidth":
      case "size":
      case "borderSize":
        option.strokeWidth = value;
        break;
      case "fontSize":
        option.fontSize = typeof value === "number" ? value + "px" : value;
      case "verticalAlign":
        if (value == "central") {
          option.dominantBaseline = "central";
        } else if (value == "bottom") {
          option.dominantBaseline = "text-bottom";
        } else if (value == "top") {
          option.dominantBaseline = "hanging";
        } else {
          option.dominantBaseline = "auto";
        }
        break;
      case "horizontalAlign":
        if (value == "right") {
          option.textAnchor = "start";
        } else if (value == "middle") {
          option.textAnchor = "middle";
        } else if (value == "end") {
          option.textAnchor = "end";
        }
        break;
      case "gradient":
      case "rotate":
      case "texture":
        break;
      case "textDirection":
        if (value == "vertical") {
          option.writingMode = "tb";
        } else if (value == "horizontal") {
          option.writingMode = "";
        }
        break;
      case "isDash":
        if (typeof value == "boolean") {
          if (value) {
            option.dasharray = "4 2";
            option.strokeDasharray = "4 2";
          } else {
            option.dashArray = "0,0";
            option.strokeDasharray = "0,0";
          }
        } else if (typeof value == "string") {
          option.dashArray = value;
        }
        default:
          option[key] = value;
    }
  }
  return option;
}

/**
 * @param {String} str 
 * @param {String} fileName 
 */
function downloadFile(str, fileName) {
  if (!fileName || "" === fileName.trim()) {
    fileName = getRandomId("IBSCharts");
  }
  let urlObj = window.URL || window.webkitURL || window;
  let url = urlObj.createObjectURL(new Blob([str]));
  downloaBaseFunc(url, fileName);
}

/**
 * 
 * @param {SVGAElement} target 
 * @param {String} fileName 
 * @param {String} type 
 */
function downloadSvg(target, fileName, type = "svg") {
  if (!fileName || "" === fileName.trim()) {
    fileName = getRandomId("IBSCharts");
  }
  let href = getSvgBase64Href(target);
  type = type.toLowerCase();
  if ("svg" === type) {
    downloaBaseFunc(href, fileName);
  } else {
    downloadPngOrJpg(href, fileName, type);
  }
}

/**
 * 
 * @param {String} href 
 * @param {String} fileName 
 */
function downloaBaseFunc(href, fileName) {
  // let aLink = document.createElementNS("http://www.w3.org/1999/xhtml", "a");
  let aLink = document.createElement("a");
  aLink.download = fileName;
  aLink.href = href;
  aLink.dispatchEvent(new MouseEvent('click', {
    bubbles: true,
    cancelable: true,
    view: window
  })); //兼容火狐
}

/**
 * 将svg转成base64字符串，并返回其href
 * @param {SVGAElement} svgElement 
 */
function getSvgBase64Href(svgElement) {
  let svgContent = svgElement.outerHTML;
  //svg转base64;
  let href = 'data:image/svg+xml;base64,' + window.btoa(unescape(encodeURIComponent(svgContent)));
  return href;
}

/**
 * 
 * @param {String} href 
 * @param {String} fileName 
 * @param {String} imgType 
 */
function downloadPngOrJpg(href, fileName, imgType) {
  let image = new Image();
  image.src = href;
  console.log("download img");
  // console.log(href);
  // console.log(image);
  image.onload = () => { // 要先确保图片完整获取到，这是个异步事件
    let canvas = document.createElement('canvas');
    console.log("load img");
    canvas.width = image.width; // 直接这样会画布的宽度和高度和原来的svg不一致，导致图片下载不全
    canvas.height = image.height;
    let context = canvas.getContext('2d');
    context.drawImage(image, 0, 0);
    let urlType = "jpg" === imgType ? "image/jpeg" : "image/png";
    // let urlType = "jpg" === imgType ? "image/tiff" : "image/png";
    // console.log(urlType);
    let imgDataUri = canvas.toDataURL(urlType, 1);
    downloaBaseFunc(imgDataUri, fileName);
  }
}


function getNow() {
  let date = new Date();
  let y = date.getFullYear();
  let m = date.getMonth() + 1;
  let d = date.getDate();
  let h = date.getHours();
  let minutes = date.getMinutes();
  return `${y}${m}${d}${h}${minutes}`;
}

/**
 * Insert an element after the targetElement
 * @param {HTMLElement} newElement 
 * @param {HTMLElement} targetElement 
 */
function insertAfter(newElement, targetElement) {
  let parent = targetElement.parentNode;
  if (parent.lastElementChild == targetElement) {
    parent.appendChild(newElement);
  } else {
    parent.insertBefore(newElement, targetElement.nextElementSibling);
  }
}

function hideAllBorder() {
  let nodeBorders = document.getElementsByClassName("nodeBorder");
  for (let i = 0; i < nodeBorders.length; i++) {
    nodeBorders[i].setAttribute("visibility", "hidden");
  }
}

function hideBorder(targetNode) {
  targetNode.style.visibility = "hidden";
}

function markerOptionConverToDomainOption(markerOption) {
  let domainOption = {
    id: markerOption.id,
    position: markerOption.position,
    shape: markerOption.shape,
    name: markerOption.text.content,
    nameTextStyle: {
      fontFamily: markerOption.text.style.fontFamily,
      fontSize: markerOption.text.style.fontSize,
      color: markerOption.text.style.color,
      rotate: markerOption.text.style.rotate,
      location: markerOption.text.position
    },
    style: {
      color: markerOption.style.color,
      texture: markerOption.style.texture,
      gradient: markerOption.style.gradient,
      rotate: markerOption.style.rotate,
      height: parseInt(markerOption.height),
    },
    borderStyle: {
      color: markerOption.style.borderColor,
      size: markerOption.style.borderSize,
      isDash: markerOption.style.isDash
    }
  }
  domainOption.position.start.site = parseInt(domainOption.position.start.site);
  domainOption.position.end.site = parseInt(domainOption.position.end.site);
  return domainOption;
}

function calculateSiteCoordinate(siteObj, parentObj) {
  let location = (siteObj.option.location + "").split(";");
  let pStartSite = parentObj.option.position.start.site;
  let pStartX = parentObj.option.coordinate.horizontal.start;
  let pStartY = parentObj.option.coordinate.vertical.start;
  let height = parentObj.option.style.height;
  let ratio = parentObj.config.lengthRatio;
  let siteCoordinate = [];
  for (let i = 0; i < location.length; i++) {
    let x = (parseFloat(location[i]) - pStartSite) * ratio + pStartX;
    let y = siteObj.option.coordinate.cy > pStartY ? pStartY + height : pStartY;
    siteCoordinate.push([x, y]);
  }
  return siteCoordinate;
}


function createShape(shape, shapeOption) {
  let content;
  switch (shape) {
    case "triangle": //
    case "rhombus": //
    case "parallelogram": //
    case "trapezium": //
    case "pentagon": //
    case "hexagon": //
    case "octagon": //
    case "arrow": //
    case "dovetailArrow": //
      content = createSvgElement("polygon", shapeOption);
      break;
      // case "circle":
      //   content = createSvgElement("ellipse", shapeOption);
      //   break;
    case "rect":
      // content = createSvgElement("rect", shapeOption);
      // break;
    case "circle":
      // content = createSvgElement("ellipse", shapeOption);
      // break;
    case "roundRect":
    case "cylinder":
      let points = shapeOption.points;
      delete shapeOption.points;
      content = createSvgElement("path", shapeOption);
      SVG(content).plot(points);
      break;
    case "none":
      content = createSvgElement("circle", shapeOption);
      break;
  }
  return content;
}

export default {
  supplementToStandardObj,
  deepCopyObj,
  getFirstElement,
  stopBubble,
  getRandomId,
  getRandomNumber,
  createSvgElement,
  getDomainLabelOption,
  createStartOrEndLabel,
  updateLabel,
  showBorderByMousedown,
  showBorderAlone,
  createBorder,
  deleteObj,
  getHexColor,
  createGradient,
  createTexture,
  downloadFile,
  downloadSvg,
  layerMove,
  Cosines,
  calAbsoluteCoordinate,
  HideAndDisplayControlBox,
  calRectEightPoints,
  calRectEightPointsWithRotate,
  standarizeStyleOption,
  getNow,
  insertAfter,
  hideAllBorder,
  hideBorder,
  markerOptionConverToDomainOption,
  calculateSiteCoordinate,
  createShape,
  calculatePointAfterRotate,
  downloaBaseFunc
}
