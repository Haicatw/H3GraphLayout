using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Matrix4d
{
    public Matrix4x4 matrix4d;

    public Matrix4d()
    {
        matrix4d = new Matrix4x4();
    }

    public void SetColumn(int index, Point4d ColumnValue)
    {
        matrix4d.SetColumn(index, ColumnValue.GetVector4());
    }

    public void SetRow(int index, Point4d ColumnValue)
    {
        matrix4d.SetRow(index, ColumnValue.GetVector4());
    }

    public Point4d GetRow(int index)
    {
        return matrix4d.GetRow(index);
    }

    public Point4d GetColumn(int index)
    {
        return matrix4d.GetColumn(index);
    }

    public void transform(Point4d v)
    {
        float x = this.matrix4d.m00 * v.x + (this.matrix4d.m01 * v.y) + (this.matrix4d.m02 * v.z) + (this.matrix4d.m03 * v.w);
        float y = this.matrix4d.m10 * v.x + (this.matrix4d.m11 * v.y) + (this.matrix4d.m12 * v.z) + (this.matrix4d.m13 * v.w);
        float z = this.matrix4d.m20 * v.x + (this.matrix4d.m21 * v.y) + (this.matrix4d.m22 * v.z) + (this.matrix4d.m23 * v.w);
        float w = this.matrix4d.m30 * v.x + (this.matrix4d.m31 * v.y) + (this.matrix4d.m32 * v.z) + (this.matrix4d.m33 * v.w);

        v.x = x; v.y = y; v.z = z; v.w = w;
    }

    public void rotX(float theta)
    {
        float cos_theta = Mathf.Cos(theta);
        float sin_theta = Mathf.Sin(theta);
        float neg_sin_theta = -sin_theta;

        this.matrix4d.m00 = 1; this.matrix4d.m01 = 0; this.matrix4d.m02 = 0; this.matrix4d.m03 = 0;
        this.matrix4d.m10 = 0; this.matrix4d.m11 = cos_theta; this.matrix4d.m12 = neg_sin_theta; this.matrix4d.m13 = 0;
        this.matrix4d.m20 = 0; this.matrix4d.m21 = sin_theta; this.matrix4d.m22 = cos_theta; this.matrix4d.m23 = 0;
        this.matrix4d.m30 = 0; this.matrix4d.m31 = 0; this.matrix4d.m32 = 0; this.matrix4d.m33 = 1;
    }

    public void rotY(float theta)
    {
        float cos_theta = Mathf.Cos(theta);
        float sin_theta = Mathf.Sin(theta);
        float neg_sin_theta = -sin_theta;

        this.matrix4d.m00 = cos_theta; this.matrix4d.m01 = 0; this.matrix4d.m02 = sin_theta; this.matrix4d.m03 = 0;
        this.matrix4d.m10 = 0; this.matrix4d.m11 = 1; this.matrix4d.m12 = 0; this.matrix4d.m13 = 0;
        this.matrix4d.m20 = neg_sin_theta; this.matrix4d.m21 = 0; this.matrix4d.m22 = cos_theta; this.matrix4d.m23 = 0;
        this.matrix4d.m30 = 0; this.matrix4d.m31 = 0; this.matrix4d.m32 = 0; this.matrix4d.m33 = 1;
    }

    public void rotZ(float theta)
    {
        float cos_theta = Mathf.Cos(theta);
        float sin_theta = Mathf.Sin(theta);
        float neg_sin_theta = -sin_theta;

        this.matrix4d.m00 = cos_theta; this.matrix4d.m01 = neg_sin_theta; this.matrix4d.m02 = 0; this.matrix4d.m03 = 0;
        this.matrix4d.m10 = sin_theta; this.matrix4d.m11 = cos_theta; this.matrix4d.m12 = 0; this.matrix4d.m13 = 0;
        this.matrix4d.m20 = 0; this.matrix4d.m21 = 0; this.matrix4d.m22 = 1; this.matrix4d.m23 = 0;
        this.matrix4d.m30 = 0; this.matrix4d.m31 = 0; this.matrix4d.m32 = 0; this.matrix4d.m33 = 1;
    }

    public void setIdentity()
    {
        this.matrix4d.m00 = 1; this.matrix4d.m01 = 0; this.matrix4d.m02 = 0; this.matrix4d.m03 = 0;
        this.matrix4d.m10 = 0; this.matrix4d.m11 = 1; this.matrix4d.m12 = 0; this.matrix4d.m13 = 0;
        this.matrix4d.m20 = 0; this.matrix4d.m21 = 0; this.matrix4d.m22 = 1; this.matrix4d.m23 = 0;
        this.matrix4d.m30 = 0; this.matrix4d.m31 = 0; this.matrix4d.m32 = 0; this.matrix4d.m33 = 1;
    }

    public static Matrix4d operator *(Matrix4d lhs, Matrix4d rhs)
    {
        return lhs.matrix4d * rhs.matrix4d;
    }

    public static Matrix4d operator *(Matrix4d lhs, float rhs)
    {
        Matrix4d temp = new Matrix4d();
        temp.SetRow(0, lhs.GetRow(0));
        temp.SetRow(1, lhs.GetRow(1));
        temp.SetRow(2, lhs.GetRow(2));
        temp.SetRow(3, lhs.GetRow(3));
        return temp;
    }

    public static implicit operator Matrix4d(Matrix4x4 m4x4)
    {
        Matrix4d temp = new Matrix4d();
        temp.SetRow(0, m4x4.GetRow(0));
        temp.SetRow(1, m4x4.GetRow(1));
        temp.SetRow(2, m4x4.GetRow(2));
        temp.SetRow(3, m4x4.GetRow(3));
        return temp;
    }
}