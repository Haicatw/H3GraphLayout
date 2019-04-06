using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GraphLayout : MonoBehaviour
{

    //public string fileName; 
    //Assets/Data/processedGraph.json
    public string networkFolder;
    public string networkName;
    public float leafRadius;
    public Node graphRoot;
    private List<Node> graphNodesList;
    //private List<GameObject> nodesPrimitives;
    //private List<GameObject> edgeHolders;
    private float EPSILON = Mathf.Epsilon;
    private Graph graphContainer;
    private Dictionary<String, int> nodeNameToIdDict;
    

    // Use this for initialization
    void Start()
    {
        //EdgeHolder.AddComponent<LineRenderer>();
        //Configration
        //leafRadius = 0.001f;

        //Initialization
        //graphRoot = new Node(); 
        graphNodesList = new List<Node>();
        //nodesPrimitives = new List<GameObject>();
        //edgeHolders = new List<GameObject>();
        nodeNameToIdDict = new Dictionary<string, int>();

        //Read in data
        ReadInGraphData();
        DebugPrint(this.graphRoot);
        CountNodeDecendents(this.graphRoot);
        //Compute initial layout
        ComputeRadius(this.graphRoot);
        SortSubtrees(this.graphRoot);
        ComputePolar(this.graphRoot);
        DebugPrint(this.graphRoot);
    }

    // Update is called once per frame
    void Update()
    {

    }

    public void readFile()
    {
        string filename = Application.streamingAssetsPath + "/" + networkFolder + "/" + networkName;
        //print(filename);
        this.graphContainer = JsonUtility.FromJson<Graph>(File.ReadAllText(filename));
    }

    public void ReadInGraphData()
    {
        readFile();
        int counter = 0;
        foreach(GraphNode node in this.graphContainer.nodes)
        {
            node.id = counter;
            this.graphNodesList.Add(new Node(node.id, node.nodeName, node.level, node.isRootNode));
            this.nodeNameToIdDict.Add(node.nodeName, node.id);

            if (node.isRootNode)
            {
                this.graphRoot = graphNodesList[counter];
            }

            counter++;
        }

        foreach(Link link in this.graphContainer.links)
        {
            //print(link.source);
            Node tempNodeA = this.graphNodesList[this.nodeNameToIdDict[link.source]];
            Node tempNodeB = this.graphNodesList[this.nodeNameToIdDict[link.target]];
            
            if (tempNodeA.nodeLevel > tempNodeB.nodeLevel)
            {
                tempNodeB.nodeChildren.Add(tempNodeA);
            } else
            {
                tempNodeA.nodeChildren.Add(tempNodeB);
            }
        }
    }

    public void CountNodeDecendents(Node parentNode)
    {
        if (parentNode.nodeChildren.Count == 0)
        {
            parentNode.nodeNumDecendents = 0;
            return;
        }
        else
        {
            int totalDecendents = 0;
            foreach(var child in parentNode.nodeChildren)
            {
                CountNodeDecendents(child);
                totalDecendents += child.nodeNumDecendents;
            }
            parentNode.nodeNumDecendents = totalDecendents;
            return;
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
            if (Math.Abs(deltaThetaCumulative) < EPSILON && Math.Abs(deltaPhiCumulatiive) < EPSILON)
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

    public void DebugPrint(Node parentNode)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            foreach(var child in parentNode.nodeChildren)
            {
                Debug.Log(child.nodeId);
                Debug.Log(child.nodeName);
                Debug.Log(child.nodeNumDecendents);
                Debug.Log(child.nodeHemspherePhi);
                Debug.Log(child.nodeHemsphereRadius);
                Debug.Log(child.nodeHemsphereTheta);

                DebugPrint(child);
                return;
            }
        }
    }

    /*

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
    */

}


public class Node
{
    public int nodeId;
    public string nodeName;
    public int nodeLevel;
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

    public Node(int nodeId, string nodeName, int nodeLevel, bool nodeIsRoot)
    {
        this.nodeId = nodeId;
        this.nodeName = nodeName;
        this.nodeLevel = nodeLevel;
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
        this.nodeIsRoot = nodeIsRoot;
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
        
        //this.nodeEuclideanPosition = HyperbolicMath.GetRotationMatrix(parentNode) * parentNode.nodeRelativeHyperbolicProjectionPosition;
        //this.nodeEuclideanPosition = HyperbolicMath.GetTranslationMatrix(parentNode) * parentNode.nodeRelativeHyperbolicProjectionPosition;
        
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

[System.Serializable]
public class Graph
{
    public Link[] links;
    public metadata graph_properties;
    public GraphNode[] nodes;
}

[System.Serializable]
public class Link
{
    public string source;
    public string target;
    public float value;
}

[System.Serializable]
public class metadata
{
    public int NumNodes;
    public int NumLinks;
    public int MaxLevel;
}

[System.Serializable]
public class GraphNode
{
    public int id;
    public string nodeName;
    public int level;
    public bool isRootNode;
}
