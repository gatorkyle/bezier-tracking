# bezier-tracking
Scalable Bezier curve defined by freely-manipulatable control points, using both mouse and keyboard input for 3-axis transformations. Created in C++ with OpenGL in Spring of 2021.

Please watch the provided .mp4 for a short demonstration on using this tool.

Keyboard Instructions:

1- Increase subdivision level of blue curve. Cycles from 0-5 before wrapping back.

2- Hide/show the red Bezier curve.

4- Toggles the double view, wherein a side view of the model can be seen below the front view. Modified from project specifications so that picking will still work in this mode.

5- Toggles the animation. In this mode, a yellow point will appear and move along the closed Bernstein-Bezier curve in both views, along with unit vectors of the point's tangent,
normal, and binormal. Note that picking is NOT enabled while the animation is running.

Shift- While the shift key is held, picking will instead move points along the Z axis. The amount of movement is determined by the mouse's horizontal motion while picking.

Mouse Instructions:

Click and drag: Move points.
