function setCam(cam)
set(gca, 'CameraPosition', cam.CameraPosition);
set(gca, 'CameraTarget', cam.CameraTarget);
set(gca, 'CameraUpVector', cam.CameraUpVector);
set(gca, 'CameraViewAngle', cam.CameraViewAngle);
