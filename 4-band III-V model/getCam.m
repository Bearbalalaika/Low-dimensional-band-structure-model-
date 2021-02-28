function cam = getCam
cam.CameraPosition = get(gca, 'CameraPosition');
cam.CameraTarget = get(gca, 'CameraTarget');
cam.CameraUpVector = get(gca, 'CameraUpVector');
cam.CameraViewAngle = get(gca, 'CameraViewAngle');

