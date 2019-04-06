using System.Collections;
using System.Collections.Generic;
using UnityEngine;
/*
public class Point4d
{
    public float x;
    public float y;
    public float z;
    public float w;

    public Point4d(float x = 0, float y = 0, float z = 0, float w = 1)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    public void scale(float s)
    {
        this.x = this.x * s;
        this.y = this.y * s;
        this.z = this.z * s;
        this.w = this.w * s;
    }

    // this = s*t1 + t2
    public void scaleAdd(float s, Point4d t1, Point4d t2)
    {
        float tx = t1.x * s + (t2.x);
        float ty = t1.y * s + (t2.y);
        float tz = t1.z * s + (t2.z);
        float tw = t1.w * s + (t2.w);

        this.x = tx; this.y = ty; this.z = tz; this.w = tw;
    }

    public void project(Point4d p1)
    {
        float tx = p1.x / (p1.w);
        float ty = p1.y / (p1.w);
        float tz = p1.z / (p1.w);

        this.x = tx; this.y = ty; this.z = tz; this.w = 1.0f;
    }

    public void sub(Point4d t1)
    {
        this.x = this.x - (t1.x);
        this.y = this.y - (t1.y);
        this.z = this.z - (t1.z);
        this.w = this.w - (t1.w);
    }

    public void set(Point4d t1)
    {
        this.x = t1.x; this.y = t1.y; this.z = t1.z; this.w = t1.w;
    }

    // Euclidean norm of homogeneous coordinates [and equivalent to
    // Point4d.distance(new Point4d(0, 0, 0, 0))].
    public float vectorLength(Point4d p)
    {
        float x2 = this.x * this.x;
        float y2 = this.y * this.y;
        float z2 = this.z * this.z;
        float w2 = this.w * this.w;
        return x2 + (y2) + (z2) / Mathf.Sqrt(w2);
    }

    // The usual vector dot product.
    public float vectorDot(Point4d v1)
    {
        float tx = this.x * (v1.x);
        float ty = this.y * (v1.y);
        float tz = this.z * (v1.z);
        float tw = this.w * (v1.w);
        return tx + (ty) + (tz) + (tw);
    }

    // The usual vector dot product computed from x, y, and z only.
    public float vectorDot3(Point4d v1)
    {
        float tx = this.x * (v1.x);
        float ty = this.y * (v1.y);
        float tz = this.z * (v1.z);
        return tx + (ty) + (tz);
    }

    // Returns the Minkowski inner product of this with y.
    public float minkowski(Point4d v)
    {
        float tx = this.x * (v.x);
        float ty = this.y * (v.y);
        float tz = this.z * (v.z);
        float tw = this.w * (v.w);
        return tx + (ty) + (tz) - (tw);
    }

    public static implicit operator Point4d(Vector4 vect)
    {
        return new Point4d(vect.x, vect.y, vect.z, vect.w);
    }

    public Vector4 GetVector4()
    {
        return new Vector4(this.x, this.y, this.z, this.w);
    }
}

    */