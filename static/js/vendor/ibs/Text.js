import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";
import IbsUtils from "./IbsUtils.js";
import { Rotater } from "./Rotater.js";

const defaultTextOption = {
  content: "",
  position: {
    x: 0,
    y: 0,
  },
  id: "",
  displayName: "",
  style: {
    fontSize: 14,
    fontFamily: "Arial",
    color: "black",
    fontStyle: "normal",
    fontWeight: "normal",
    textDecoration: "none", //textDecoration
    horizontalAlign: "middle", // textAnchor
    verticalAlign: "center", //[center, buttom, top]
    rotate: 0,
    textDirection: "horizontal", //
  },
  borderStyle: {
    //外部矩形边框
    color: "none", //stroke
    size: 1, //strokeWidth
    isDash: false, //dashArray=1,1
  },
  editable: true,
};
/**
 *
 * @param {{content:String, id:String, displayName: String, position:{x:Number,y:Number}, style:{fontSize:Number, fontFamily:String,
 * color:String, textStyle:String, textWeight:String, textDecoration:String, horizontalAlign:String, verticalAlign: String, rotate: Number},
 * borderStyle:{color:String, size: Number, isDash:Boolean},editable: Boolean}} TextOption
 * @param {{IbsCharts}} parent
 */
function Text(TextOption, parent) {
  this.parent = parent;
  let { option, config } = this.standardizeConfigAndOption(TextOption);
  this.option = option; //设定默认值
  this.config = config;
  this.children = [];
  this.nodesIndex = null;
  if (parent && "root" === parent.config.type) {
    this.nodesIndex = parent.children.length;
  }
  this.nodes = this.createNodes();
  // if (this.option.editable) {
  // let nbox = SVG(this.config.textNode).bbox();
  // let rotaterOption = {
  //   x: parseFloat(nbox.cx) - 6,
  //   y: parseFloat(nbox.y) - 20,
  //   rotate: this.option.style.rotate
  // }
  // if (this.config.parentType != "marker") {
  //   var rotater = new Rotater(rotaterOption, this);
  //   this.children[0] = rotater;
  //   this.nodes.appendChild(rotater.nodes);
  //   this.nodes.appendChild(rotater.config.cloneNode);
  //  //IbsUtils.HideAndDisplayControlBox(rotater.nodes, this.nodes);
  // }
  // }
  this.updateTransform();
  this.draggable(this.option.editable);
}
/**
 * standarize the input text option
 * @param {{}} textOption
 */
Text.prototype.standardizeConfigAndOption = function (
  textOption,
  updateId = true
) {
  let option = IbsUtils.supplementToStandardObj(defaultTextOption, textOption);
  if (updateId) {
    option.id = IbsUtils.getRandomId("Text");
  }
  option.displayName = option.displayName ? option.displayName : option.id;
  option.style.fontSize = Math.max(0, option.style.fontSize);
  option.borderStyle.size = Math.max(0, option.borderStyle.size);
  let config = {};
  config.nodesTransform = [0, 0, 0, option.position.x, option.position.y];
  let parentType = this.parent.config.type;
  config.parentType = parentType;
  config.type = "text";
  return {
    option,
    config,
  };
};
/**
 * create nodes of text object
 */
Text.prototype.createNodes = function () {
  let text = IbsUtils.createSvgElement("text", {
    x: this.option.position.x,
    y: this.option.position.y,
    id: this.option.id,
  });
  let styleObj = IbsUtils.standarizeStyleOption(this.option.style);
  for (let key in styleObj) {
    if ("stroke" !== key) {
      text.style[key] = styleObj[key];
    }
  }
  text.textContent = this.option.content;
  text.style.userSelect = "none";
  this.config.textNode = text;
  let group = IbsUtils.createSvgElement("g", {
    id: this.option.id + "Group",
  });
  group.appendChild(text);
  if (this.option.borderStyle.color != "none") {
    let border = this.createBorder();
    group.appendChild(border);
  }
  if (this.option.editable) {
    let container = this.createContainer();
    this.config.containerNode = container;
    group.appendChild(container);
    group.appendChild(this.children[0].config.cloneNode); //将rotater的cloneNode添加到group中
    IbsUtils.showBorderByMousedown(group, container);
    IbsUtils.showBorderAlone(container);
  }
  switch (this.parent.config.type) {
    case "marker":
    case "site":
    case "domain":
      break;
    case "root":
      this.parent.config.selected = this;
      group.addEventListener("mousedown", () => {
        this.parent.config.isHideAll = false;
        this.parent.config.selected = this;
      });
  }
  text.removeAttribute("svgjs:data");
  return group;
};
/**
 * create border node
 */
Text.prototype.createBorder = function () {
  let nbox = SVG(this.config.textNode).bbox();
  let borderAttr = {
    x: nbox.x - 5,
    y: nbox.y - 2,
    width: nbox.width + 10,
    height: nbox.height + 4,
    style: `fill:none;stroke:${this.option.borderStyle.color};stroke-width:${this.option.borderStyle.size};`,
  };
  let border = IbsUtils.createSvgElement("rect", borderAttr);
  if (this.option.borderStyle.isDash) {
    border.style.strokeDasharray = "5,5";
  }
  this.config.borderNode = border;
  return border;
};
/**
 * create the rotater node
 */
Text.prototype.createContainer = function () {
  let nbox = SVG(this.config.textNode).bbox();
  let { x, y, width, height } = nbox;
  let container = IbsUtils.createSvgElement("g", { class: "nodeBorder" });
  for (let i = 0; i < 2; i++) {
    for (let j = 0; j < 2; j++) {
      let cx = parseFloat(x) + parseInt(j) * parseFloat(width);
      let cy = parseFloat(y) + parseInt(i) * parseFloat(height);
      let pointAttr = {
        cx,
        cy,
        r: 4,
        fill: "#ffa500",
      };
      let point = IbsUtils.createSvgElement("circle", pointAttr);
      container.appendChild(point);
    }
  }
  let rotaterOption = {
    x: this.option.position.x - 6,
    y: this.option.position.y - 30,
    rotate: this.option.style.rotate,
  };
  let rotater = new Rotater(rotaterOption, this);
  this.children.push(rotater);
  container.appendChild(rotater.nodes);
  container.setAttribute("visibility", "hidden");
  return container;
};
/**
 * make the node draggable
 * @param {Boolean} editable
 */
Text.prototype.draggable = function (editable = true) {
  if (editable) {
    var textSVG = SVG(this.nodes);
    textSVG.draggable();
    var dy; //text对应的box的cy与文字基线的差值
    var dx;
    var targetX;
    var targetY;
    var nbox;
    var self = this;
    textSVG.on(`beforedrag.Text${self.option.id}`, (e) => {
      if ("root" === this.parent.config.type) {
        this.parent.config.isUndo = false;
        let operation = {
          target: this,
          cmd: "drag",
          args: JSON.stringify(this.option),
          nodesIndex: this.nodesIndex,
        };
        this.parent.config.undoStack.push(operation);
      }
      nbox = SVG(self.nodes).bbox();
      dy = self.option.position.y - nbox.cy; //text对应的box的cy与文字基线的差值
      dx = self.option.position.x - nbox.cx;
    });
    textSVG.on(`dragmove.Text${self.option.id}`, (e) => {
      let { box } = e.detail;
      e.preventDefault();
      targetX = Math.round((parseFloat(box.cx) + dx) * 100) / 100;
      targetY = Math.round((parseFloat(box.cy) + dy) * 100) / 100;
      self.Move(targetX, targetY);
      //self.Move(targetX, targetY);
    });
    textSVG.on(`dragend.Text${self.option.id}`, () => {
      if (self.children.length) {
        self.children[0].updateCloneNode();
      }
    });
  }
};

/**
 * update the text node based on the text option
 * @param {{}} TextOption
 */
Text.prototype.update = function (TextOption) {
  let previousOption = JSON.stringify(this.option);
  let { option, config } = this.standardizeConfigAndOption(TextOption);
  this.option = option;
  this.config = config;
  this.children = [];
  let updateNode = this.createNodes();
  this.parent.nodes.replaceChild(updateNode, this.nodes);
  this.nodes = updateNode;
  // if (TextOption.editable) {
  //   let nbox = SVG(this.config.textNode).bbox();
  //   let rotaterOption = {
  //     x: nbox.cx - 6,
  //     y: nbox.y - 20,
  //     rotate: this.option.style.rotate
  //   }
  //   if (this.config.parentType != "Markers") {
  //     var rotater = new Rotater(rotaterOption, this);
  //     this.children[0] = rotater;
  //   }
  //   this.nodes.appendChild(rotater.nodes);
  //   this.nodes.appendChild(rotater.config.cloneNode);
  // }
  this.draggable(this.option.editable);
  this.updateTransform();
  // IbsUtils.HideAndDisplayControlBox(this.children[0].nodes,this.nodes,this.parent.nodes);
  // IbsUtils.HideAndDisplayControlBox(this.children[0].config.cloneNode,this.nodes,this.parent.nodes);
  if ("root" === this.parent.config.type) {
    this.parent.config.isUndo = false;
    let operation = {
      target: this,
      cmd: "update",
      args: previousOption,
      nodesIndex: this.nodesIndex,
    };
    this.parent.config.undoStack.push(operation);
  }
  return this;
};
Text.prototype.delete = function () {
  IbsUtils.deleteObj(this);
  if ("root" === this.parent.config.type) {
    // this.parent.config.elementsCount["text"] -= 1;
    this.parent.config.isUndo = false;
    let operation = {
      target: this,
      cmd: "delete",
      args: JSON.stringify(this.option),
      nodesIndex: this.nodesIndex,
    };
    this.parent.config.undoStack.push(operation);
  }
  return true;
};
/**
 * update the innerHTML of the node
 * @param {String} newContent new innerHTML
 * @param {Boolean} adjust  change this.content or not
 */
Text.prototype.updateContent = function (newContent, adjust = true) {
  this.config.textNode.textContent = newContent;
};
/**
 * move the group to a new position
 * @param {Number} input_x
 * @param {Number} input_y
 */
Text.prototype.Move = function (input_x, input_y, notSite = true) {
  let translateX =
    typeof input_x == "number"
      ? (translateX = input_x - this.option.position.x)
      : false;
  let translateY =
    typeof input_y == "number"
      ? (translateY = input_y - this.option.position.y)
      : false;
  this.option.position.x = input_x;
  this.option.position.y = input_y;
  if (notSite) {
    this.updateTransform(translateX, translateY, null);
  }
};
/**
 * update the transform attribute of nodes
 * @param {Number} x
 * @param {Number} y
 * @param {Number} angel
 * @param {Number} cx
 * @param {Number} cy
 */
Text.prototype.updateTransform = function (x, y, angel, cx, cy) {
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
  if (typeof cx == "number") {
    this.config.nodesTransform[3] = cx;
  }
  if (typeof cy == "number") {
    this.config.nodesTransform[4] = cy;
  }
  var transformValue = `translate(${this.config.nodesTransform[0]},${this.config.nodesTransform[1]}),rotate(${this.option.style.rotate},${this.config.nodesTransform[3]},${this.config.nodesTransform[4]})`;
  this.config.textNode.setAttribute("transform", transformValue);
  // for (let i in this.children) {
  //   this.children[i].nodes.setAttribute("transform", transformValue);
  // }
  if (this.option.borderStyle.color != "none") {
    this.config.borderNode.setAttribute("transform", transformValue);
  }
  if (this.config.containerNode) {
    this.config.containerNode.setAttribute("transform", transformValue);
  }
};
export { Text };
