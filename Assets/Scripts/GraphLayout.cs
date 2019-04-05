using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Graph : MonoBehaviour
{

    [SerializeField] string fileName; //Assets/Data/processedGraph.json
    [SerializeField] float leafRadius;
    public Node graphRoot;
    public List<Node> graphNodesList;
    public List<GameObject> nodesPrimitives;
    public List<GameObject> edgeHolders;

    // Use this for initialization
    void Start()
    {
        //EdgeHolder.AddComponent<LineRenderer>();
        //Configration
        leafRadius = 0.001f;

        //Initialization
        //graphRoot = new Node(); 
        graphNodesList = new List<Node>();
        nodesPrimitives = new List<GameObject>();
        edgeHolders = new List<GameObject>();

        //Read in data
        ReadInGraphData("");
        //Compute initial layout
        ComputeRadius(graphRoot);
        SortSubtrees(graphRoot);
        ComputePolar(graphRoot);

    }

    private void ReadInGraphData(string path)
    {
        if (System.IO.File.Exists(path))
        {

        }
    }

    public void SortSubtrees(Node parentNode)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            parentNode.nodeChildren.Sort((x, y) => y.nodeNumDecendents.CompareTo(x.nodeNumDecendents));
            foreach (var child in parentNode.nodeChildren)
            {
                SortSubtrees(child);
            }
        }
        return;
    }

    // Update is called once per frame
    void Update()
    {

    }

    public void ComputeRadius(Node parentNode)
    {
        //bottom up algorithm. (From leaf to root)
        //With given leaf size, calculate parent hemsphere size and recursively compute upper level hemsphere size
        float parentHemsphereArea = 0.0f; //Hp = sum of pi*r^2
        if (parentNode.nodeNumDecendents == 0)
        {
            parentNode.nodeHemsphereRadius = leafRadius;
            return;
        }
        else
        {
            foreach (var child in parentNode.nodeChildren)
            {
                ComputeRadius(child); //recursive bottom up call
                //Euclidean Space Case for understnding:
                //parentHemsphereArea += Mathf.Pow(child.nodeHemsphereRadius, 2);
                //Hyperbolic Space Case:
                parentHemsphereArea += (float)(Math.Cosh(child.nodeHemsphereRadius) - 1);
            }
        }
        //Euclidean Space Case:
        //parentNode.nodeHemsphereRadius = Mathf.Sqrt(parentHemsphereArea);
        //Hyperbolic Space Case:
        parentNode.nodeHemsphereRadius = (float)Math.Sinh(parentHemsphereArea);
        return;
    }

    public void ComputePolar(Node parentNode)
    {
        //TODO: hyperbolic space translation and rotation by their parent
        float deltaThetaCumulative = 0.0f;
        float deltaPhiCumulatiive = 0.0f;
        float currentGreatestChildRadius = 0.0f;

        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }

        foreach (var child in parentNode.nodeChildren)
        {
            if (deltaThetaCumulative == 0 && deltaPhiCumulatiive == 0)
            {
                child.nodeHemspherePhi = 0.0f;
                child.nodeHemsphereTheta = 0.0f;
                currentGreatestChildRadius = parentNode.nodeHemsphereRadius;
            }
            else
            {
                if (deltaThetaCumulative > (Mathf.PI * 2))
                {
                    deltaPhiCumulatiive += (float)Math.Atan(Math.Tanh(currentGreatestChildRadius) / Math.Sinh(parentNode.nodeHemsphereRadius));
                    deltaThetaCumulative = 0.0f;
                    currentGreatestChildRadius = child.nodeHemsphereRadius;//children is already sorted
                    //TODO: check if first hemisphere need exception
                }

                child.nodeHemspherePhi = deltaPhiCumulatiive;
                child.nodeHemsphereTheta = (float)Math.Atan(Math.Tanh(child.nodeHemsphereRadius) / (Math.Sinh(parentNode.nodeHemsphereRadius) * Math.Sin(deltaPhiCumulatiive)));
            }
            ComputePolar(child);
        }
    }

    public void ComputeGlobalPolar(Node parentNode)
    {
        foreach (var child in parentNode.nodeChildren)
        {
            if (child.nodeNumDecendents == 0)
            {
                child.CalculateGlobalCartesianCoordinate(parentNode);
                return;
            }
            else
            {
                child.CalculateGlobalCartesianCoordinate(parentNode);
                ComputeGlobalPolar(parentNode);
            } 
        }
        return;
    }

    public void DrawNodes()
    {
        nodesPrimitives.Add(GameObject.CreatePrimitive(PrimitiveType.Sphere));
        nodesPrimitives[nodesPrimitives.Count - 1].transform.position = new Vector4(graphRoot.nodeEuclideanPosition.x, graphRoot.nodeEuclideanPosition.y, graphRoot.nodeEuclideanPosition.z, graphRoot.nodeEuclideanPosition.w);
        DrawNodesRecursive(graphRoot);
    }

    public void DrawNodesRecursive(Node node)
    {
        foreach (var child in node.nodeChildren)
        {
            nodesPrimitives.Add(GameObject.CreatePrimitive(PrimitiveType.Sphere));
            nodesPrimitives[nodesPrimitives.Count - 1].transform.position = new Vector4(graphRoot.nodeEuclideanPosition.x, graphRoot.nodeEuclideanPosition.y, graphRoot.nodeEuclideanPosition.z, graphRoot.nodeEuclideanPosition.w);
            DrawNodesRecursive(child);
        }
    }

    public void DrawEdges()
    {
        //TODO: deprecated, not sure GC rules
        edgeHolders.Add(new GameObject());
        edgeHolders[edgeHolders.Count - 1].AddComponent<LineRenderer>();
        LineRenderer lr = edgeHolders[edgeHolders.Count - 1].GetComponent<LineRenderer>();
        lr.material = new Material(Shader.Find("Particles/Alpha Blended Premultiply"));
        //lr.SetColors(color, color);
        //TODO: set colors
        lr.SetWidth(0.1f, 0.1f);
        lr.SetPosition(0, new Vector4(graphRoot.nodeEuclideanPosition.x, graphRoot.nodeEuclideanPosition.y, graphRoot.nodeEuclideanPosition.z, graphRoot.nodeEuclideanPosition.w));
        foreach (var child in graphRoot.nodeChildren)
        {
            //lr.SetPosition(1, );
        }
    }

    public void DrawEdgesRecursive(Node node)
    {

    }

}


public class Node
{
    public Point4d nodeEuclideanPosition;
    public Point4d nodeRelativeHyperbolicProjectionPosition;
    public float nodeAnglePhi;
    public float nodeAngleTheta;
    public float nodeR;
    public float nodeHemspherePhi;
    public float nodeHemsphereTheta;
    public float nodeHemsphereRadius;
    public List<Node> nodeChildren;
    public Node nodeParent;
    public int nodeNumDecendents;
    public bool nodeIsRoot;

    public Node()
    {
        this.nodeEuclideanPosition = new Point4d(0, 0, 0, 1);
        this.nodeRelativeHyperbolicProjectionPosition = new Point4d(0, 0, 0, 1);
        this.nodeAnglePhi = 0.0f;
        this.nodeAngleTheta = 0.0f;
        this.nodeR = 0.0f;
        this.nodeHemspherePhi = 0.0f;
        this.nodeHemsphereTheta = 0.0f;
        this.nodeHemsphereRadius = 0.0f;
        this.nodeChildren = new List<Node>();
        this.nodeParent = null;
        this.nodeNumDecendents = 0;
        this.nodeIsRoot = false;
    }

    public void SetNodeEuclideanPosition(float x, float y, float z, float w)
    {
        this.nodeEuclideanPosition = new Point4d(x, y, z, w);
    }

    public void SetNodeRelativeHyperbolicProjectionPosition(float x, float y, float z, float w)
    {
        this.nodeRelativeHyperbolicProjectionPosition = new Point4d(x, y, z, w);
    }

    public void CalculateAndSetNodeEuclideanPosition()
    {
        float x = this.nodeR * Mathf.Sin(this.nodeAnglePhi) * Mathf.Cos(this.nodeAngleTheta);
        float y = this.nodeR * Mathf.Sin(this.nodeAnglePhi) * Mathf.Sin(this.nodeAngleTheta);
        float z = this.nodeR * Mathf.Cos(this.nodeAnglePhi);

        this.SetNodeEuclideanPosition(x, y, z, 1);
    }

    public void CalculateAndSetNodeRelativeHyperbolicProjectionPosition()
    {
        float x = this.nodeHemsphereRadius * Mathf.Sin(this.nodeHemspherePhi) * Mathf.Cos(this.nodeHemsphereTheta);
        float y = this.nodeHemsphereRadius * Mathf.Sin(this.nodeHemspherePhi) * Mathf.Sin(this.nodeHemsphereTheta);
        float z = this.nodeHemsphereRadius * Mathf.Cos(this.nodeHemspherePhi);

        this.SetNodeRelativeHyperbolicProjectionPosition(x, y, z, 1);
    }

    public void SetHemspheres(float Phi, float Theta)
    {
        this.nodeHemspherePhi = Phi;
        this.nodeHemsphereTheta = Theta;
    }

    public void SetPolarCoordinates(float R, float Phi, float Theta)
    {
        this.nodeR = R;
        this.nodeAnglePhi = Phi;
        this.nodeAngleTheta = Theta;
    }

    internal void CalculateGlobalCartesianCoordinate(Node parentNode)
    {
        
        this.nodeEuclideanPosition = HyperbolicMath.GetRotationMatrix(parentNode) * parentNode.nodeRelativeHyperbolicProjectionPosition;
        this.nodeEuclideanPosition = HyperbolicMath.GetTranslationMatrix(parentNode) * parentNode.nodeRelativeHyperbolicProjectionPosition;
        
    }
}

class HyperbolicMath
{

    public static float MinkowskiInnerProduct(Point4d x, Point4d y)
    {
        return (x.x * y.x + x.y * y.y + x.z * y.z - x.w * y.w);
    }

    public static Matrix4d GetReflectionMatrix(Point4d point)
    {
        //rflection = I - 2 * p * pTrans * I31 / <p, p>h
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

public class TransformUtility
{

}
