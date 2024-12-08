#version 330 core
precision mediump float;

out vec4 fragColor;

void main() {
    fragColor = vec4(gl_FragCoord.xy / vec2(800, 600), 0.0, 1.0);
}