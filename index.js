/**
 * This is Barnes Hut simulation algorithm for 3d case. Implementation
 * is highly optimized (avoids recusion and gc pressure)
 *
 * http://www.cs.princeton.edu/courses/archive/fall03/cs126/assignments/barnes-hut.html
 *
 * NOTE: This module duplicates a lot of code from 2d case. Primary reason for
 * this is performance. Every time I tried to abstract away vector opertaions
 * I had negative impact on performance. So in this case I'm sacrifying code
 * reuse in favor of speed
 */

module.exports = function(options) {
  options = options || {};
  options.gravity = typeof options.gravity === 'number' ? options.gravity : -1;
  options.theta = typeof options.theta === 'number' ? options.theta : 0.8;

  // we require deterministic randomness here
  var random = require('ngraph.random').random(1984),
    Node = require('./node'),
    InsertStack = require('./insertStack'),
    isSamePosition = require('./isSamePosition');

  var gravity = options.gravity,
    updateQueue = [],
    insertStack = new InsertStack(),
    theta = options.theta,

    nodesCache = [],
    currentInCache = 0,
    newNode = function() {
      // To avoid pressure on GC we reuse nodes.
      var node = nodesCache[currentInCache];
      if (node) {
        node.quads[0] = null;
        node.quads[4] = null;
        node.quads[1] = null;
        node.quads[5] = null;
        node.quads[2] = null;
        node.quads[6] = null;
        node.quads[3] = null;
        node.quads[7] = null;
        node.body = null;
        node.mass = node.massX = node.massY = node.massZ = 0;
        node.left = node.right = node.top = node.bottom = node.front = node.back = 0;
      } else {
        node = new Node();
        nodesCache[currentInCache] = node;
      }

      ++currentInCache;
      return node;
    },

    root = newNode(),

    // Inserts body to the tree
    insert = function(newBody) {
      insertStack.reset();
      insertStack.push(root, newBody);

      while (!insertStack.isEmpty()) {
        var stackItem = insertStack.pop(),
          node = stackItem.node,
          body = stackItem.body;

        if (!node.body) {
          // This is internal node. Update the total mass of the node and center-of-mass.
          var x = body.pos.x;
          var y = body.pos.y;
          var z = body.pos.z;
          node.mass += body.mass;
          node.massX += body.mass * x;
          node.massY += body.mass * y;
          node.massZ += body.mass * z;

          // Recursively insert the body in the appropriate quadrant.
          // But first find the appropriate quadrant.
          var quadIdx = 0, // Assume we are in the 0's quad.
            left = node.left,
            right = (node.right + left) / 2,
            top = node.top,
            bottom = (node.bottom + top) / 2,
            back = node.back,
            front = (node.front + back) / 2;

          if (x > right) { // somewhere in the eastern part.
            quadIdx += 1;
            var oldLeft = left;
            left = right;
            right = right + (right - oldLeft);
          }
          if (y > bottom) { // and in south.
            quadIdx += 2;
            var oldTop = top;
            top = bottom;
            bottom = bottom + (bottom - oldTop);
          }
          if (z > front) { // and in frontal part
            quadIdx += 4;
            var oldBack = back;
            back = front;
            front = back + (back - oldBack);
          }

          var child = node.quads[quadIdx];
          if (!child) {
            // The node is internal but this quadrant is not taken. Add subnode to it.
            child = newNode();
            child.left = left;
            child.top = top;
            child.right = right;
            child.bottom = bottom;
            child.back = back;
            child.front = front;
            child.body = body;

            node.quads[quadIdx] = child;
          } else {
            // continue searching in this quadrant.
            insertStack.push(child, body);
          }
        } else {
          // We are trying to add to the leaf node.
          // We have to convert current leaf into internal node
          // and continue adding two nodes.
          var oldBody = node.body;
          node.body = null; // internal nodes do not cary bodies

          if (isSamePosition(oldBody.pos, body.pos)) {
            // Prevent infinite subdivision by bumping one node
            // anywhere in this quadrant
            var retriesCount = 3;
            do {
              var offset = random.nextDouble();
              var dx = (node.right - node.left) * offset;
              var dy = (node.bottom - node.top) * offset;
              var dz = (node.front - node.back) * offset;

              oldBody.pos.x = node.left + dx;
              oldBody.pos.y = node.top + dy;
              oldBody.pos.z = node.back + dz;
              retriesCount -= 1;
              // Make sure we don't bump it out of the box. If we do, next iteration should fix it
            } while (retriesCount > 0 && isSamePosition(oldBody.pos, body.pos));

            if (retriesCount === 0 && isSamePosition(oldBody.pos, body.pos)) {
              // This is very bad, we ran out of precision.
              // if we do not return from the method we'll get into
              // infinite loop here. So we sacrifice correctness of layout, and keep the app running
              // Next layout iteration should get larger bounding box in the first step and fix this
              return;
            }
          }
          // Next iteration should subdivide node further.
          insertStack.push(node, oldBody);
          insertStack.push(node, body);
        }
      }
    },

    update = function(sourceBody) {
      var queue = updateQueue,
        v,
        dx, dy, dz,
        r, fx = 0,
        fy = 0,
        fz = 0,
        queueLength = 1,
        shiftIdx = 0,
        pushIdx = 1;

      queue[0] = root;

      while (queueLength) {
        var node = queue[shiftIdx],
          body = node.body;

        queueLength -= 1;
        shiftIdx += 1;
        // technically there should be external "if (body !== sourceBody) {"
        // but in practice it gives slightghly worse performance, and does not
        // have impact on layout correctness
        if (body && body !== sourceBody) {
          // If the current node is a leaf node (and it is not source body),
          // calculate the force exerted by the current node on body, and add this
          // amount to body's net force.
          dx = body.pos.x - sourceBody.pos.x;
          dy = body.pos.y - sourceBody.pos.y;
          dz = body.pos.z - sourceBody.pos.z;
          r = Math.sqrt(dx * dx + dy * dy + dz * dz);

          if (r === 0) {
            // Poor man's protection against zero distance.
            dx = (random.nextDouble() - 0.5) / 50;
            dy = (random.nextDouble() - 0.5) / 50;
            dz = (random.nextDouble() - 0.5) / 50;
            r = Math.sqrt(dx * dx + dy * dy + dz * dz);
          }

          // This is standard gravition force calculation but we divide
          // by r^3 to save two operations when normalizing force vector.
          v = gravity * body.mass * sourceBody.mass / (r * r * r);
          fx += v * dx;
          fy += v * dy;
          fz += v * dz;
        } else {
          // Otherwise, calculate the ratio s / r,  where s is the width of the region
          // represented by the internal node, and r is the distance between the body
          // and the node's center-of-mass
          dx = node.massX / node.mass - sourceBody.pos.x;
          dy = node.massY / node.mass - sourceBody.pos.y;
          dz = node.massZ / node.mass - sourceBody.pos.z;
          r = Math.sqrt(dx * dx + dy * dy + dz * dz);

          if (r === 0) {
            // Sorry about code duplucation. I don't want to create many functions
            // right away. Just want to see performance first.
            dx = (random.nextDouble() - 0.5) / 50;
            dy = (random.nextDouble() - 0.5) / 50;
            dz = (random.nextDouble() - 0.5) / 50;
            r = Math.sqrt(dx * dx + dy * dy + dz * dz);
          }
          // If s / r < θ, treat this internal node as a single body, and calculate the
          // force it exerts on sourceBody, and add this amount to sourceBody's net force.
          if ((node.right - node.left) / r < theta) {
            // in the if statement above we consider node's width only
            // because the region was squarified during tree creation.
            // Thus there is no difference between using width or height.
            v = gravity * node.mass * sourceBody.mass / (r * r * r);
            fx += v * dx;
            fy += v * dy;
            fz += v * dz;
          } else {
            // Otherwise, run the procedure recursively on each of the current node's children.

            // I intentionally unfolded this loop, to save several CPU cycles.
            if (node.quads[0]) {
              queue[pushIdx] = node.quads[0];
              queueLength += 1;
              pushIdx += 1;
            }
            if (node.quads[1]) {
              queue[pushIdx] = node.quads[1];
              queueLength += 1;
              pushIdx += 1;
            }
            if (node.quads[2]) {
              queue[pushIdx] = node.quads[2];
              queueLength += 1;
              pushIdx += 1;
            }
            if (node.quads[3]) {
              queue[pushIdx] = node.quads[3];
              queueLength += 1;
              pushIdx += 1;
            }
            if (node.quads[4]) {
              queue[pushIdx] = node.quads[4];
              queueLength += 1;
              pushIdx += 1;
            }
            if (node.quads[5]) {
              queue[pushIdx] = node.quads[5];
              queueLength += 1;
              pushIdx += 1;
            }
            if (node.quads[6]) {
              queue[pushIdx] = node.quads[6];
              queueLength += 1;
              pushIdx += 1;
            }
            if (node.quads[7]) {
              queue[pushIdx] = node.quads[7];
              queueLength += 1;
              pushIdx += 1;
            }
          }
        }
      }

      sourceBody.force.x += fx;
      sourceBody.force.y += fy;
      sourceBody.force.z += fz;
    },

    insertBodies = function(bodies) {
      var x1 = Number.MAX_VALUE,
        y1 = Number.MAX_VALUE,
        z1 = Number.MAX_VALUE,
        x2 = Number.MIN_VALUE,
        y2 = Number.MIN_VALUE,
        z2 = Number.MIN_VALUE,
        i,
        max = bodies.length;

      // To reduce quad tree depth we are looking for exact bounding box of all particles.
      i = max;
      while (i--) {
        var pos = bodies[i].pos;
        var x = pos.x;
        var y = pos.y;
        var z = pos.z;
        if (x < x1) {
          x1 = x;
        }
        if (x > x2) {
          x2 = x;
        }
        if (y < y1) {
          y1 = y;
        }
        if (y > y2) {
          y2 = y;
        }
        if (z < z1) {
          z1 = z;
        }
        if (z > z2) {
          z2 = z;
        }
      }

      // Squarify the bounds.
      var maxSide = Math.max(x2 - x1, Math.max(y2 - y1, z2 - z1));

      x2 = x1 + maxSide;
      y2 = y1 + maxSide;
      z2 = z1 + maxSide;

      currentInCache = 0;
      root = newNode();
      root.left = x1;
      root.right = x2;
      root.top = y1;
      root.bottom = y2;
      root.back = z1;
      root.front = z2;

      i = max - 1;
      if (i > 0) {
        root.body = bodies[i];
      }
      while (i--) {
        insert(bodies[i], root);
      }
    };

  return {
    insertBodies: insertBodies,
    updateBodyForce: update,
    options: function(newOptions) {
      if (newOptions) {
        if (typeof newOptions.gravity === 'number') {
          gravity = newOptions.gravity;
        }
        if (typeof newOptions.theta === 'number') {
          theta = newOptions.theta;
        }

        return this;
      }

      return {
        gravity: gravity,
        theta: theta
      };
    }
  };
};
