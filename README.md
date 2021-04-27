# bezier-tracking
Scalable Bezier curve defined by freely-manipulatable control points, using both mouse and keyboard input for 3-axis transformations. Created in C++ with OpenGL in Spring of 2021.

![Demo](https://user-images.githubusercontent.com/42983161/116212935-56f44b80-a713-11eb-9c36-bdea91ad3d3c.gif)

    Keyboard Instructions:
  
    1- Increase subdivision level of blue curve. Cycles from 0-5 before wrapping back.
  
    2- Hide/show the red Bezier curve.

    4- Toggles the double view, wherein a side view of the model can be seen below the front view. Modified from project specifications so that picking will still work in this mode.

    5- Toggles the animation. In this mode, a yellow point will appear and move along the closed Bernstein-Bezier curve in both views, along with unit vectors of the point's tangent, normal, and binormal. Note that picking is NOT enabled while the animation is running.

    Shift- While the shift key is held, picking will instead move control points along the Z axis. The amount of movement is determined by the mouse's horizontal motion while picking.

    Mouse Instructions:

    Click and drag: Move the control points.
