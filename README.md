# Pies For Maya

Pies for Maya is a plugin integration of our soft-body physics engine Pies into Autodesk Maya - checkout the [Pies repository](https://github.com/nithinp7/Pies) to learn more about its features and implementation. This plugin introduces easy-to-use tools to create, constrain, and simulate soft-body physics meshes. 

## Highlights

<img src="https://github.com/nithinp7/PiesForMaya/blob/main/Media/CoilingRope.gif">
<img src="https://github.com/nithinp7/PiesForMaya/blob/main/Media/JellyGears.gif">
<img src="https://github.com/nithinp7/PiesForMaya/blob/main/Media/CoilGears120.gif">

##### A simple scene with a number of objects falling past a position-constrained object.
<img src="https://github.com/nithinp7/PiesForMaya/blob/main/Media/BasicScene.gif">

##### In the below video, a simple keyframed motion of a Pies Constraint drags along a soft-body.
<img src="https://github.com/nithinp7/PiesForMaya/blob/main/Media/KeyframedConstraints.gif">

##### In the below video, a keyframed Pies Constraint pulls the torus upwards, crashing into other falling objects.
<img src="https://github.com/nithinp7/PiesForMaya/blob/main/Media/KeyframedConstraints2.gif">

## Instructions

The `Pies` menu on the top-bar allows users to quickly create a solver, convert meshes into soft-bodies, and constrain / link objects. The behavior of the solver and the soft-bodies can be further fine-tuned in the attribute editor:
- Adjust the solver iteration count, collision error, gravity, damping, friction, and much more on the `PiesSolverNode`.
- Adjust stiffness parameters, strain limits, volume preservance, density, and initial velocity for each soft-body on the corresponding `PiesSoftBodyNode`
- Click `Add Link` in the `Pies` menu to place shape-matching links in the scene that rigidly connect all particles in the specified radius.
- Click `Add Constraint` in the `Pies` menu to place constraints that fix all particles in the specified radius at the specified position and orientation. This constraint can be keyframed and move / rotate over time, as shown in previous videos.
