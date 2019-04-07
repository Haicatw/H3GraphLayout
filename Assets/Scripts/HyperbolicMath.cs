using System.Collections;
using System.Collections.Generic;
using UnityEngine;
/*
class HyperbolicMath
{

    public static float MinkowskiInnerProduct(Point4d x, Point4d y)
    {
        return (x.x * y.x + x.y * y.y + x.z * y.z - x.w * y.w);
    }

    public static Matrix4d GetReflectionMatrix(Point4d point)
    {
        //rflection = I - 2 * p * pTrans * I31 / <p, p>h
        
        float xx = point.x * point.x;
        float xy = point.x * point.y;
        float xz = point.x * point.z;
        float xw = point.x * point.w;

        float yy = point.y * point.y;
        float yz = point.y * point.z;
        float yw = point.y * point.w;

        float zz = point.z * point.z;
        float zw = point.z * point.w;

        float ww = point.w * point.w;

        float pp_h = xx + yy + zz - ww;
        float temp = -2.0f / pp_h;

        Matrix4d ppTI31 = new Matrix4d();

        ppTI31.SetRow(0, new Point4d(xx * temp + 1, xy * temp, xz * temp, -xw * temp));
        ppTI31.SetRow(1, new Point4d(xy * temp, yy * temp + 1, yz * temp, -yw * temp));
        ppTI31.SetRow(2, new Point4d(xz * temp, yz * temp, zz * temp + 1, -zw * temp));
        ppTI31.SetRow(3, new Point4d(xw * temp, yw * temp, zw * temp, -ww * temp + 1));

        return ppTI31;
    }

    public static Point4d GetMidpoint(Point4d pointA, Point4d pointB)
    {
        float coefficientA = Mathf.Sqrt(MinkowskiInnerProduct(pointB, pointB) * MinkowskiInnerProduct(pointA, pointB));
        float coefficientB = Mathf.Sqrt(MinkowskiInnerProduct(pointA, pointA) * MinkowskiInnerProduct(pointA, pointB));

        return new Point4d(pointA.x * coefficientA + pointB.x * coefficientB, pointA.y * coefficientA + pointB.y * coefficientB, pointA.z * coefficientA + pointB.z * coefficientB, pointA.w * coefficientA + pointB.w * coefficientB);
    }

    public static float dotProduct(Vector3 x, Vector3 y)
    {
        return x.x * y.x + x.y * y.y + x.z * y.z;
    }

    public static float dotProduct(Point4d x, Point4d y)
    {
        return (x.x * y.x + x.y * y.y + x.z * y.z) / (x.w * y.w);
    }

    public static float vectorLength(Point4d p)
    {
        float w2 = p.w * p.w;
        return Mathf.Sqrt((p.x * p.x + p.y * p.y + p.z * p.z) / w2);
    }

    public static void makeUnitVector(Point4d p)
    {
        float s = p.w * vectorLength(p);
        p.x /= s;
        p.y /= s;
        p.z /= s;
        p.w = 1.0f;
    }


}
*/

/*
    Matrix4x4 I = new Matrix4x4();
    I.SetRow(0, new Vector4(1, 0, 0, 0));
    I.SetRow(1, new Vector4(0, 1, 0, 0));
    I.SetRow(2, new Vector4(0, 0, 1, 0));
    I.SetRow(3, new Vector4(0, 0, 0, 1));

    Matrix4x4 I31 = new Matrix4x4();
    I31.SetRow(0, new Vector4(1, 0, 0, 0));
    I31.SetRow(1, new Vector4(0, 1, 0, 0));
    I31.SetRow(2, new Vector4(0, 0, 1, 0));
    I31.SetRow(3, new Vector4(0, 0, 0, -1));

    Matrix4x4 PointMultiplyPointTranspose = new Matrix4x4();
    //PointMultiplyPointTranspose.SetRow(0, new Vector4(point.x * point.x * 2, point.x * point.y * 2, point.x * point.z * 2, point.x * point.w * 2));
    //PointMultiplyPointTranspose.SetRow(1, new Vector4(point.x * point.y * 2, point.y * point.y * 2, point.y * point.z * 2, point.y * point.w * 2));
    //PointMultiplyPointTranspose.SetRow(2, new Vector4(point.x * point.z * 2, point.y * point.z * 2, point.z * point.z * 2, point.z * point.w * 2));
    //PointMultiplyPointTranspose.SetRow(3, new Vector4(point.x * point.w * 2, point.y * point.w * 2, point.z * point.w * 2, point.w * point.w * 2));
    //TODO: cross product

    float minkowski = MinkowskiInnerProduct(point, point);

    Matrix4x4 result = new Matrix4x4();

    //TODO: check if it is cross product or dot product
    result = PointMultiplyPointTranspose * I31;

    result.SetRow(0, new Vector4(result[0, 0] / minkowski, result[0, 1] / minkowski, result[0, 2] / minkowski, result[0, 3] / minkowski));
    result.SetRow(1, new Vector4(result[1, 0] / minkowski, result[1, 1] / minkowski, result[1, 2] / minkowski, result[1, 3] / minkowski));
    result.SetRow(2, new Vector4(result[2, 0] / minkowski, result[2, 1] / minkowski, result[2, 2] / minkowski, result[2, 3] / minkowski));
    result.SetRow(3, new Vector4(result[3, 0] / minkowski, result[3, 1] / minkowski, result[3, 2] / minkowski, result[3, 3] / minkowski));

    result.SetRow(0, new Vector4(1 - result[0, 0], 0 - result[0, 1], 0 - result[0, 2], 0 - result[0, 3]));
    result.SetRow(1, new Vector4(0 - result[1, 0], 1 - result[1, 1], 0 - result[1, 2], 0 - result[1, 3]));
    result.SetRow(2, new Vector4(0 - result[2, 0], 0 - result[2, 1], 1 - result[2, 2], 0 - result[2, 3]));
    result.SetRow(3, new Vector4(0 - result[3, 0], 0 - result[3, 1], 0 - result[3, 2], 1 - result[3, 3]));

    return result;
    */
