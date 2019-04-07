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
    private List<GameObject> nodesPrimitives;
    private List<GameObject> edgeHolders;
    private float EPSILON = Mathf.Epsilon;
    private Graph graphContainer;
    private Dictionary<String, int> nodeNameToIdDict;
    private float globalScaler;
    

    // Use this for initialization
    void Start()
    {
        //EdgeHolder.AddComponent<LineRenderer>();
        //Configration
        //leafRadius = 0.001f;

        //Initialization
        //graphRoot = new Node(); 
        graphNodesList = new List<Node>();
        nodesPrimitives = new List<GameObject>();
        edgeHolders = new List<GameObject>();
        nodeNameToIdDict = new Dictionary<string, int>();
        globalScaler = 10.0f;

        //Read in data
        ReadInGraphData();
        //DebugPrint(this.graphRoot);
        Debug.Log("Finish Loading");
        CountNodeDecendents(this.graphRoot);
        Debug.Log("Finish Counting");
        //Compute initial layout
        ComputeRadius(this.graphRoot);
        Debug.Log("Finish Calculating Radius");
        SortSubtrees(this.graphRoot);
        Debug.Log("Finish Sorting");
        ComputePolar(this.graphRoot);
        ComputeRelativeHyperbolicProjectionPosition(this.graphRoot);
        Debug.Log("Finish Polar");
        ComputeCoordinatesEuclidean(this.graphRoot);
        //ComputeCoordinates();
        DrawNodes();
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
                tempNodeA.nodeParent = tempNodeB;
            } else
            {
                tempNodeA.nodeChildren.Add(tempNodeB);
                tempNodeB.nodeParent = tempNodeA;
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
            foreach(Node child in parentNode.nodeChildren)
            {
                CountNodeDecendents(child);
                totalDecendents += child.nodeNumDecendents;
                totalDecendents += 1; //count the child itself
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
            foreach (Node child in parentNode.nodeChildren)
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
        
        if (parentNode.nodeNumDecendents == 0)
        {
            parentNode.nodeHemsphereRadius = leafRadius;
            return;
        }
        else
        {
            float parentHemsphereArea = 0.0f; //Hp = sum of pi*r^2
            foreach (Node child in parentNode.nodeChildren)
            {
                ComputeRadius(child); //recursive bottom up call
                //Euclidean Space Case for understnding:
                //parentHemsphereArea += Mathf.Pow(child.nodeHemsphereRadius, 2);
                //Hyperbolic Space Case:
                parentHemsphereArea += (float)(Math.Cosh(child.nodeHemsphereRadius) - 1);
            }
            //Euclidean Space Case:
            //parentNode.nodeHemsphereRadius = Mathf.Sqrt(parentHemsphereArea);
            //Hyperbolic Space Case:
            parentNode.nodeHemsphereRadius = (float)HyperbolicMath.ASinh(Math.Sqrt(parentHemsphereArea));
        }
        
        return;
    }

    public void ComputePolar(Node parentNode)
    {
        float deltaThetaCumulative = 0.0f;
        float deltaPhiCumulatiive = 0.0f;
        float currentGreatestChildRadius = 0.0f;

        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }

        foreach (Node child in parentNode.nodeChildren)
        {
            if (Math.Abs(deltaThetaCumulative) < EPSILON && Math.Abs(deltaPhiCumulatiive) < EPSILON)
            {
                child.nodeHemspherePhi = 0.0f;
                child.nodeHemsphereTheta = 0.0f;
                currentGreatestChildRadius = parentNode.nodeHemsphereRadius;
                deltaThetaCumulative = 2 * Mathf.PI * 2; //any constant that greater than 2pi would work here
            }
            else
            {
                deltaThetaCumulative += (float)Math.Atan(Math.Tanh(child.nodeHemsphereRadius) / (Math.Sinh(parentNode.nodeHemsphereRadius) * Math.Sin(deltaPhiCumulatiive)));

                if (deltaThetaCumulative > (Mathf.PI * 2))
                {
                    deltaPhiCumulatiive += (float)Math.Atan(Math.Tanh(currentGreatestChildRadius) / Math.Sinh(parentNode.nodeHemsphereRadius));
                    deltaThetaCumulative = 0.0f;
                    currentGreatestChildRadius = child.nodeHemsphereRadius;//children is already sorted
                    //TODO: check if first hemisphere need exception
                } 
                
                child.nodeHemspherePhi = deltaPhiCumulatiive;
                child.nodeHemsphereTheta = deltaThetaCumulative;
            }
            ComputePolar(child);
        }
    }

    public void ComputeGlobalPolar(Node parentNode)
    {
        //TODO: hyperbolic space translation and rotation by their parent
        foreach (Node child in parentNode.nodeChildren)
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

    public void ComputeRelativeHyperbolicProjectionPosition(Node parentNode)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            foreach(Node child in parentNode.nodeChildren)
            {
                child.CalculateAndSetNodeRelativeHyperbolicProjectionPosition();
                ComputeRelativeHyperbolicProjectionPosition(child);
            }
        }
    }

    public void DebugPrint(Node parentNode)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            foreach(Node child in parentNode.nodeChildren)
            {
                Debug.Log("Id: " + child.nodeId);
                Debug.Log("Name: " + child.nodeName);
                Debug.Log("Coord: " + child.nodeEuclideanPosition.x + " " + child.nodeEuclideanPosition.y + " " + child.nodeEuclideanPosition.z );

                DebugPrint(child);
            }
            return;
        }
    }

    private void ComputeCoordinates()
    {
        // The root node is always positioned at the origin.
        this.graphRoot.nodeEuclideanPosition = HyperbolicTransformation.ORIGIN4;
        ComputeCoordinatesSubtree(HyperbolicTransformation.I4, graphRoot);
    }

    public void ComputeCoordinatesSubtree(Matrix4d parentTransform, Node parentNode)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }

        float parentRadiusE = HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius);

        float lastPhi = 0.0f;
        Matrix4d rotPhi = HyperbolicTransformation.I4;

        Point4d childCenterAbsolute = new Point4d();
        Point4d childPoleAbsolute = new Point4d();
        foreach (Node child in parentNode.nodeChildren)
        {
            float childRadiusE = HyperbolicMath.euclideanDistance(child.nodeHemsphereRadius);
            float childPhi = child.nodeHemspherePhi;

            if (!(Mathf.Abs(childPhi - lastPhi) < EPSILON)) 
            {
                lastPhi = childPhi;
                rotPhi = HyperbolicTransformation.buildZRotation(childPhi);
            }

            Matrix4d rot = HyperbolicTransformation.buildXRotation(child.nodeHemsphereTheta);

            rot = rot * rotPhi;

            childCenterAbsolute.set(parentRadiusE, 0.0f, 0.0f, 1.0f);
            rot.transform(childCenterAbsolute);
            float childPoleE = HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius + child.nodeHemsphereRadius);

            childPoleAbsolute.set(childPoleE, 0.0f, 0.0f, 1.0f);
            rot.transform(childPoleAbsolute);

            parentTransform.transform(childCenterAbsolute);
            parentTransform.transform(childPoleAbsolute);

            child.nodeEuclideanPosition = childCenterAbsolute;
            //graph.setNodeLayoutCoordinates(child, childCenterAbsolute);

            Matrix4d childTransform = HyperbolicTransformation.buildCanonicalOrientation(childCenterAbsolute, childPoleAbsolute);

            ComputeCoordinatesSubtree(childTransform, child);
        }
    }

    public void ComputeCoordinatesEuclidean(Node parentNode)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {

            foreach(Node child in parentNode.nodeChildren)
            {
                Quaternion rotation = Quaternion.Euler(child.nodeHemsphereTheta, 0.0f, child.nodeHemspherePhi);
                // Set the translation, rotation and scale parameters.
                Matrix4x4 m = Matrix4x4.TRS(new Vector3(HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius), 0.0f, 0.0f), rotation, new Vector3(globalScaler, globalScaler, globalScaler));
                Vector3 tempVect3 = m.MultiplyPoint3x4(Point4d.projectAndGetVect3(child.nodeRelativeHyperbolicProjectionPosition));
                child.nodeEuclideanPosition = new Point4d(tempVect3.x, tempVect3.y, tempVect3.z);
                ComputeCoordinatesEuclidean(child);
            }
        }
    }

    public void DrawNodes()
    {
        nodesPrimitives.Add(GameObject.CreatePrimitive(PrimitiveType.Sphere));
        nodesPrimitives[nodesPrimitives.Count - 1].transform.position = new Vector3(graphRoot.nodeEuclideanPosition.x * 10000, graphRoot.nodeEuclideanPosition.y * 10000, graphRoot.nodeEuclideanPosition.z * 10000);
        nodesPrimitives[nodesPrimitives.Count - 1].name = graphRoot.nodeName;
        DrawNodesRecursive(graphRoot);
    }

    public void DrawNodesRecursive(Node node)
    {
        foreach (Node child in node.nodeChildren)
        {
            nodesPrimitives.Add(GameObject.CreatePrimitive(PrimitiveType.Sphere));
            nodesPrimitives[nodesPrimitives.Count - 1].transform.position = new Vector3(child.nodeEuclideanPosition.x * 10000, child.nodeEuclideanPosition.y * 10000, child.nodeEuclideanPosition.z * 10000);
            nodesPrimitives[nodesPrimitives.Count - 1].name = child.nodeName;
            DrawNodesRecursive(child);
        }
    }
    /*
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
        foreach (Node child in graphRoot.nodeChildren)
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

    public Matrix4d(float m00 = 0.0f, float m01 = 0.0f, float m02 = 0.0f, float m03 = 0.0f,
              float m10 = 0.0f, float m11 = 0.0f, float m12 = 0.0f, float m13 = 0.0f,
              float m20 = 0.0f, float m21 = 0.0f, float m22 = 0.0f, float m23 = 0.0f,
              float m30 = 0.0f, float m31 = 0.0f, float m32 = 0.0f, float m33 = 0.0f)
    {
        this.matrix4d = new Matrix4x4();
        this.matrix4d.m00 = m00; this.matrix4d.m01 = m01; this.matrix4d.m02 = m02; this.matrix4d.m03 = m03;
        this.matrix4d.m10 = m10; this.matrix4d.m11 = m11; this.matrix4d.m12 = m12; this.matrix4d.m13 = m13;
        this.matrix4d.m20 = m20; this.matrix4d.m21 = m21; this.matrix4d.m22 = m22; this.matrix4d.m23 = m23;
        this.matrix4d.m30 = m30; this.matrix4d.m31 = m31; this.matrix4d.m32 = m32; this.matrix4d.m33 = m33;
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

    public void project()
    {
        float tx = this.x / (this.w);
        float ty = this.y / (this.w);
        float tz = this.z / (this.w);

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

    public void set(float x, float y, float z, float w)
    {
        this.x = x; this.y = y; this.z = z; this.w = w;
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

    public static Vector3 projectAndGetVect3(Point4d p)
    {
        float tx = p.x / (p.w);
        float ty = p.y / (p.w);
        float tz = p.z / (p.w);
        return new Vector3(tx, ty, tz);
    }
}

class HyperbolicMath
{

    public static float MinkowskiInnerProduct(Point4d x, Point4d y)
    {
        return (x.x * y.x + x.y * y.y + x.z * y.z - x.w * y.w);
    }

    public static float euclideanDistance(float x)
    {
        float y = (float)(Math.Cosh(x / 2.0));
        return Mathf.Sqrt(1.0f - 1.0f / (y * y));
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

    public static double ASinh(double value)
    {
        double temp = value + Math.Sqrt(Math.Pow(value, 2.0) + 1);
        return Math.Log(temp, Math.E);
    }
}

public class HyperbolicTransformation
{
    public static Point4d ORIGIN4 = new Point4d( 0.0f, 0.0f, 0.0f, 1.0f );
    public static Matrix4d I4 = new Matrix4d(1.0f, 0.0f, 0.0f, 0.0f,
                                 0.0f, 1.0f, 0.0f, 0.0f,
                                 0.0f, 0.0f, 1.0f, 0.0f,
                                 0.0f, 0.0f, 0.0f, 1.0f);

    public static Matrix4d buildCanonicalOrientation(Point4d a, Point4d b)
    {
        /* local scratch variables; will be transformed */
        Point4d pa = new Point4d(a.x, a.y, a.z, a.w);
        Point4d pb = new Point4d(b.x, b.y, b.z, b.w);

        Point4d pivot = findPivotPoint(pa, pb);

        Matrix4d retval = buildTranslation(ORIGIN4, pivot);

        Matrix4d t1 = buildTranslation(pivot, ORIGIN4);
        t1.transform(pa);
        t1.transform(pb);

        retval *= buildTranslation(ORIGIN4, pa);

        Matrix4d t2 = buildTranslation(pa, ORIGIN4);
        t2.transform(pa);
        t2.transform(pb);

        /* calculate spherical coordinates (rho, phi, theta) of pb */

        // Projection to affine coordinates is necessary so that we can
        // directly reference the x, y, and z components in the following
        // calculations.
        pb.project();

        float rho = HyperbolicMath.vectorLength(pb);
        float phi = Mathf.Acos(pb.x / rho);
        float theta = Mathf.Atan2(pb.z, pb.y);

        if (phi == 0.0f)
        {
            /* rotate line to achieve alignment on positive x-axis */
            retval *= buildXRotation(theta);
            retval *= buildZRotation(phi);
        }

        return retval;
    }

    public static Matrix4d buildXRotation(float angle)
    {
        Matrix4d m = new Matrix4d();
        m.rotX(angle);
        return m;
    }

    public static Matrix4d buildYRotation(float angle)
    {
        Matrix4d m = new Matrix4d();
        m.rotY(angle);
        return m;
    }

    public static Matrix4d buildZRotation(float angle)
    {
        Matrix4d m = new Matrix4d();
        m.rotZ(angle);
        return m;
    }

    public static Matrix4d buildTranslation(Point4d source, Point4d dest)
    {
        float aa_h = HyperbolicMath.MinkowskiInnerProduct(source, source);
        float bb_h = HyperbolicMath.MinkowskiInnerProduct(dest, dest);
        float ab_h = HyperbolicMath.MinkowskiInnerProduct(source, dest);
        float sourceScale = Mathf.Sqrt(bb_h * ab_h);
        float destScale = Mathf.Sqrt(aa_h * ab_h);
        Point4d midpoint = new Point4d();
        midpoint.x = sourceScale * source.x + destScale * dest.x;
        midpoint.y = sourceScale * source.y + destScale * dest.y;
        midpoint.z = sourceScale * source.z + destScale * dest.z;
        midpoint.w = sourceScale * source.w + destScale * dest.w;

        Matrix4d r_a = buildReflection(source);
        Matrix4d r_m = buildReflection(midpoint);
        r_m *= r_a;
        return r_m;
    }

    public static Matrix4d buildReflection(Point4d point)
    {
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

    private static Point4d findPivotPoint(Point4d a4, Point4d b4)
    {
        Vector3 a = new Vector3();
        Vector3 b = new Vector3();

        a = Point4d.projectAndGetVect3(a4);
        b = Point4d.projectAndGetVect3(b4);

        Vector3 a_minus_b = new Vector3(a.x - b.x, a.y - b.y, a.z - b.z);

        float p = HyperbolicMath.dotProduct(a, a_minus_b);
        float q = HyperbolicMath.dotProduct(b, a_minus_b);
        float r = HyperbolicMath.dotProduct(a_minus_b, a_minus_b);

        return new Point4d(p * b.x - q * a.x,
                   p * b.y - q * a.y,
                   p * b.z - q * a.z,
                   r);
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
