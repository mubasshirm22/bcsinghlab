import IbsUtils from "./IbsUtils.js";
import { SVG } from "./svg.js/main.js";
import "./svg.js/svg.draggable.js";

const defaultPoint = {
  cx: "",
  cy: "",
  color: "black",
  r: 3,
  index: "",
  editPoint: true
};
/**
 * create a points object
 * @param {{cx:Number, cy:Number, color:String, r:Number, index:Number, editPoint:Boolean}} pointOption
 * @param {IbsCharts} parent
 */
function Points(pointOption, parent) {
  this.option = IbsUtils.supplementToStandardObj(defaultPoint, pointOption);
  this.parent = parent;
  this.nodes = this.createNodes();
  this.draggable(this.option.editPoint);
}

/**
 * create point node based on the input option
 */
Points.prototype.createNodes = function() {
  var pointAttr = {
    cx: this.option.cx,
    cy: this.option.cy,
    r: this.option.r,
    id: this.parent.option.id + "Point" + this.option.index,
    // style: `fill:${this.option.color};stroke:none;stroke-width:10;stroke-opacity:0;`,
    style: `fill:${this.option.color};stroke:#000000;stroke-width:10;stroke-opacity:0;`
  };
  if (this.parent.config.type == "container") {
    pointAttr.style = `fill:${this.option.color};`;
  }
  var point = IbsUtils.createSvgElement("circle", pointAttr);
  return point;
};
/**
 * make the node draggable
 * @param {Boolean} editPoint
 */
Points.prototype.draggable = function(editPoint) {
  if (editPoint) {
    var self = this;
    let points = SVG(self.nodes);
    let domainMin;
    let domainMax;
    let nucleotide;
    let ref;
    let parentType = self.parent.config.type;
    points.draggable();
    points.on(`beforedrag.namespace`, e => {
      if (this.parent.parent.config.type == "domain") {
        ref = SVG(this.parent.parent.config.ref).rbox(SVG(this.nodes));
        nucleotide = SVG(
          "#" + this.parent.parent.parent.config.groupId + "_rect"
        ).rbox(SVG(this.nodes));
        // console.log(nucleotide);
        domainMin = nucleotide.x;
        domainMax = domainMin + nucleotide.width;
      }

      let target;
      let IbsRoot;
      if ("container" === parentType) {
        // marker, site, domain
        if ("marker" === self.parent.parent.config.type) {
          target = self.parent.parent;
          IbsRoot = target.parent;
        } else {
          // site, domain
          target = self.parent.parent;
          IbsRoot = target.parent.parent;
        }
      } else {
        // line, bracket
        target = self.parent;
        IbsRoot = target.parent;
      }
      IbsRoot.config.isUndo = false;

      let operation = {
        target: target,
        cmd: "reshape",
        args: JSON.stringify(target.option)
        // nodesIndex
      };
      // console.log(IbsRoot);
      IbsRoot.config.undoStack.push(operation);
      // console.log("reshape operation");
      // console.log(IbsRoot.config.undoStack);
    });
    points.on(
      `dragmove.${self.parent.option.id + "Point" + self.option.index}`,
      e => {
        let { handler, box } = e.detail;
        e.preventDefault();
        var x;
        var y;
        // let parentType = self.parent.config.type
        if (parentType == "line") {
          self.parent.reshape(self.option.index, box.cx, box.cy, true);
          handler.move(box.x, box.y);
        } else if (parentType == "container") {
          switch (self.option.index) {
            case 1:
              x = self.option.cx;
              y =
                box.cy > self.parent.config.centerPoint[1]
                  ? self.parent.config.centerPoint[1]
                  : box.cy;
              break;
            case 5:
              x = self.option.cx;
              y =
                box.cy < self.parent.config.centerPoint[1]
                  ? self.parent.config.centerPoint[1]
                  : box.cy;
              break;
            case 3:
            case 7:
              x = box.cx;
              y = self.option.cy;
              break;
            case 0:
            case 2:
              x = box.cx;
              y =
                box.cy > self.parent.config.centerPoint[1]
                  ? self.parent.config.centerPoint[1]
                  : box.cy;
              break;
            case 4:
            case 6:
              x = box.cx;
              y =
                box.cy < self.parent.config.centerPoint[1]
                  ? self.parent.config.centerPoint[1]
                  : box.cy;
              break;
          }
          //禁止拖出边界
          // if (this.parent.parent.config.type == "domain") {
          //   if (domainMin > x) {
          //     x = domainMin;
          //   } else if (domainMax < x) {
          //     x = domainMax;
          //   }
          // }
          self.parent.reshape(self.option.index, x, y);
          handler.move(x - self.option.r, y - self.option.r);
        } else if (parentType == "bracket") {
          let x;
          switch (parseInt(this.option.index)) {
            case 0:
              if (box.x > this.parent.config.nodeCoor[0]) {
                x = this.parent.config.nodeCoor[0];
              } else {
                x = box.cx;
              }
              this.parent.reshape(parseInt(this.option.index), x);
              break;
            case 1:
              if (box.x < this.parent.config.nodeCoor[0]) {
                x = this.parent.config.nodeCoor[0];
              } else {
                x = box.cx;
              }
              this.parent.reshape(
                parseInt(this.option.index),
                2 * this.parent.config.nodeCoor[0] - x
              );
              break;
          }
          handler.move(x - 3);
        }
      }
    );
    points.on(
      `dragend.${self.parent.option.id + "Point" + self.option.index}`,
      e => {
        self.updateCoordinate(true);
      }
    );
  }
};
/**
 * update Node attribute
 * @param {Number} index the point node index
 * @param {Number} cx the new X coordinate
 * @param {Number} cy the new Y coordinate
 */
Points.prototype.updateCoordinate = function(selfPoint = true, cx, cy) {
  if (selfPoint) {
    let nbox = SVG(this.nodes).bbox();
    this.option.cx = nbox.cx;
    this.option.cy = nbox.cy;
  } else {
    if (cx) {
      this.nodes.setAttribute("cx", cx);
    }
    if (cy) {
      this.nodes.setAttribute("cy", cy);
    }
  }
};

export { Points };
