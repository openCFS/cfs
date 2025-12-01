# Test cases for aligning the heat conduction tensor with a vector field

The rotation is validated for 2D and 3D cases.

The vector fields and the rotation results can be streamed to a .txt file.\
The logging streams should look like this:

```
<logging>
  <stream name="coeffunctionRotationValidation" level="dbg3">
    <file name="coeffunctionRot.txt"/>
  </stream >
</logging>
```

This .txt file can then be read by the script *validation.py*.

There the same rotation is calculated by a scipy method and the relative error of the rotation results is printed.
