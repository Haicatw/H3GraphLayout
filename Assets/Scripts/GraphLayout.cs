using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GraphLayout : MonoBehaviour
{

    [SerializeField] string fileName; //Assets/Data/processedGraph.json
    [SerializeField] float leafRadius;
    public Node graphRoot;
    public List<Node> graphNodesList;
    public List<GameObject> nodesPrimitives;
    public List<GameObject> edgeHolders;
    public float EPSILON = Mathf.Epsilon;

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
        
        //this.nodeEuclideanPosition = HyperbolicMath.GetRotationMatrix(parentNode) * parentNode.nodeRelativeHyperbolicProjectionPosition;
        //this.nodeEuclideanPosition = HyperbolicMath.GetTranslationMatrix(parentNode) * parentNode.nodeRelativeHyperbolicProjectionPosition;
        
    }
}

public class TransformUtility
{

}
