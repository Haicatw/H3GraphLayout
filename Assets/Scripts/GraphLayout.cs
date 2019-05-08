﻿using System;
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
    public double leafRadius;
    public Node graphRoot;
    public Color color;
    public Material lineMaterial;
    public float globalScaler;
    public float lineWidth;
    public float randeringRadius;
    private List<Node> graphNodesList;
    private List<GameObject> nodesPrimitives;
    private List<GameObject> edgeHolders;
    private readonly double EPSILON = Mathf.Epsilon;
    private Graph graphContainer;
    private Dictionary<String, int> nodeNameToIdDict;
    //private double globalScaler;
    

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
        //globalScaler = 1000.0;

        //Point4d a = new Point4d(0.2, 0, 0, 1);
        Point4d a = new Point4d(0, 0, 0, 1);
        Point4d b = new Point4d(-0.5, -0.5, 0, 1);
        Point4d bp = new Point4d(0.3, -0.7, 0, 1);
        Matrix4d translationMat = new Matrix4d();
        translationMat = HyperbolicMath.getTranslationMatrix(b, bp);
        Point4d at = new Point4d(0.2, 0, 0, 1);
        translationMat.transform(at);

        Point4d coord = new Point4d(0, 0, 0, 1);
        Point4d origin = new Point4d(0, 0, 0, 1);
        Matrix4d transformMat = HyperbolicMath.getTranslationMatrix(origin, new Point4d(0.5, 0, 0, 1));
        Matrix4d rotationMat = HyperbolicMath.getRotationMatrix(origin, new Point4d(0, 1, 0, 1), Math.PI / 3);
        rotationMat *= HyperbolicMath.getRotationMatrix(origin, new Point4d(0, 0, 1, 1), Math.PI / 6);
        transformMat.transform(coord);
        rotationMat.transform(coord);
        //GameObject testOb = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        //testOb.transform.position = new Vector3((float)coord.x,
        //                                        (float)coord.y,
        //                                        (float)coord.z);


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
        //Scaling(graphRoot);
        DrawNodes();
        DrawEdges(this.graphRoot);
        //DebugPrint(this.graphRoot);
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
            double parentHemsphereArea = 0.0; //Hp = sum of pi*r^2
            foreach (Node child in parentNode.nodeChildren)
            {
                ComputeRadius(child); //recursive bottom up call
                //Euclidean Space Case for understnding:
                //parentHemsphereArea += Math.Pow(child.nodeHemsphereRadius, 2);
                //Hyperbolic Space Case:
                //H3Math.TWO_PI * (H3Math.cosh(r / K) - 1.0);
                //parentHemsphereArea += (double)(Math.PI * 2 * (Math.Cosh(child.nodeHemsphereRadius / 2) - 1));
                parentHemsphereArea += (double)(Math.PI * 2 * (Math.Cosh(child.nodeHemsphereRadius) - 1));
            }
            //Euclidean Space Case:
            //parentNode.nodeHemsphereRadius = Math.Sqrt(parentHemsphereArea);
            //Hyperbolic Space Case:
            //K * H3Math.asinh(Math.sqrt(area / (H3Math.TWO_PI * K * K)));
            //parentNode.nodeHemsphereRadius = (double)HyperbolicMath.ASinh(Math.Sqrt(parentHemsphereArea / (Math.PI * 2 * 4))) * 2;
            parentNode.nodeHemsphereRadius = (double)HyperbolicMath.ASinh(Math.Sqrt(parentHemsphereArea / (Math.PI * 2)));
        }
        
        return;
    }

    public void ComputePolar(Node parentNode)
    {
        //double deltaTheta = 0.0;
        //double deltaPhi = 0.0;
        double deltaThetaCumulative = 0.0;
        double deltaPhiCumulatiive = 0.0;
        double currentGreatestChildRadius = 0.0;

        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }

        foreach (Node child in parentNode.nodeChildren)
        {
            if (Math.Abs(deltaThetaCumulative) < EPSILON && Math.Abs(deltaPhiCumulatiive) < EPSILON)
            {
                child.nodeHemspherePhi = 0.0;
                child.nodeHemsphereTheta = 0.0;
                currentGreatestChildRadius = parentNode.nodeHemsphereRadius;
                deltaThetaCumulative = 2 * Math.PI * 2; //any constant that greater than 2pi would work here
            }
            else
            {
                //Math.atan(H3Math.tanh(rn / K) / (H3Math.sinh(rp / K) * Math.sin(phi)));
                //deltaTheta = (double)Math.Atan(Math.Tanh(child.nodeHemsphereRadius / 2) / (Math.Sinh(parentNode.nodeHemsphereRadius / 2) * Math.Sin(deltaPhiCumulatiive)));
                deltaThetaCumulative += (double)Math.Atan(Math.Tanh(child.nodeHemsphereRadius) / (Math.Sinh(parentNode.nodeHemsphereRadius) * Math.Sin(deltaPhiCumulatiive)));
                //if (deltaThetaCumulative + 2*deltaTheta > (Math.PI * 2))
                if (deltaThetaCumulative > (Math.PI * 2))
                {
                    //Math.atan(H3Math.tanh(rj / K) / H3Math.sinh(rp / K));
                    //deltaPhi = (double)Math.Atan(Math.Tanh(currentGreatestChildRadius / 2) / Math.Sinh(parentNode.nodeHemsphereRadius / 2));
                    deltaPhiCumulatiive += (double)Math.Atan(Math.Tanh(currentGreatestChildRadius) / Math.Sinh(parentNode.nodeHemsphereRadius));
                    deltaThetaCumulative = 0.0;
                    currentGreatestChildRadius = child.nodeHemsphereRadius;//children is already sorted
                    //TODO: check if first hemisphere need exception
                }

                child.nodeHemspherePhi = deltaPhiCumulatiive;
                child.nodeHemsphereTheta = deltaThetaCumulative;
                
                //child.nodeHemspherePhi = deltaPhiCumulatiive + deltaPhi;
                //child.nodeHemsphereTheta = deltaThetaCumulative + deltaTheta;

                //deltaPhiCumulatiive += (2 * deltaPhi);
                //deltaThetaCumulative += (2 * deltaTheta);
            }
            ComputePolar(child);
        }
    }

    public void ComputeGlobalPolarAngle(Node parentNode)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            foreach(Node child in parentNode.nodeChildren)
            {
                //child.nodeAnglePhi += parentNode.
            }
        }
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
                child.CalculateAndSetNodeRelativeHyperbolicProjectionPosition(parentNode.nodeHemsphereRadius);
                
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
    /*
    public void ComputeCoordinates()
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

        double parentRadiusE = HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius);

        double lastPhi = 0.0;
        Matrix4d rotPhi = HyperbolicTransformation.I4;

        Point4d childCenterAbsolute = new Point4d();
        Point4d childPoleAbsolute = new Point4d();
        foreach (Node child in parentNode.nodeChildren)
        {
            double childRadiusE = HyperbolicMath.euclideanDistance(child.nodeHemsphereRadius);
            double childPhi = child.nodeHemspherePhi;

            if (!(Math.Abs(childPhi - lastPhi) < EPSILON)) 
            {
                lastPhi = childPhi;
                rotPhi = HyperbolicTransformation.buildZRotation(childPhi);
            }

            Matrix4d rot = HyperbolicTransformation.buildXRotation(child.nodeHemsphereTheta);

            rot = rot * rotPhi;

            childCenterAbsolute.set(parentRadiusE, 0.0, 0.0, 1.0);
            rot.transform(childCenterAbsolute);
            double childPoleE = HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius + child.nodeHemsphereRadius);

            childPoleAbsolute.set(childPoleE, 0.0, 0.0, 1.0);
            rot.transform(childPoleAbsolute);

            parentTransform.transform(childCenterAbsolute);
            parentTransform.transform(childPoleAbsolute);

            child.nodeEuclideanPosition = childCenterAbsolute;
            //graph.setNodeLayoutCoordinates(child, childCenterAbsolute);

            Matrix4d childTransform = new Matrix4d();//HyperbolicTransformation.buildCanonicalOrientation(childCenterAbsolute, childPoleAbsolute);

            ComputeCoordinatesSubtree(childTransform, child);
        }
    }
    */
    public void ComputeCoordinatesEuclidean(Node parentNode)
    {
        /*
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            Matrix4x4 translationMatrix = Matrix4x4.Translate(new Vector3(HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius) + parentNode.nodeEuclideanPosition.x, parentNode.nodeEuclideanPosition.y, parentNode.nodeEuclideanPosition.z));
            foreach (Node child in parentNode.nodeChildren)
            {
                Vector3 childCoord = Point4d.projectAndGetVect3(child.nodeRelativeHyperbolicProjectionPosition);
                childCoord = translationMatrix.MultiplyPoint(childCoord);

                Quaternion rotation = Quaternion.Euler(child.nodeHemsphereTheta, 0.0, child.nodeHemspherePhi);
                // Set the translation, rotation and scale parameters.
                //Matrix4x4 m = Matrix4x4.TRS(new Vector3(HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius) + parentNode.nodeEuclideanPosition.x, parentNode.nodeEuclideanPosition.y, parentNode.nodeEuclideanPosition.z), rotation, new Vector3(globalScaler, globalScaler, globalScaler));
                Matrix4x4 m = Matrix4x4.TRS(new Vector3(0.0, 0.0, 0.0), rotation, new Vector3(globalScaler, globalScaler, globalScaler));
                Vector3 tempVect3 = m.MultiplyPoint3x4(childCoord);
                child.nodeEuclideanPosition = new Point4d(tempVect3.x, tempVect3.y, tempVect3.z);
                ComputeCoordinatesEuclidean(child);
            }
        }
        */
        /*
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            //Matrix4d rotationOfCoordinate = HyperbolicTransformation.buildCanonicalOrientationEuclidean(HyperbolicTransformation.ORIGIN4, parentNode.nodeEuclideanPosition);
            double parentPhi = Point4d.getPhiByPoint(parentNode.nodeEuclideanPosition);
            double parentTheta = Point4d.getThetaByPoint(parentNode.nodeEuclideanPosition);
            Debug.Log(parentNode.nodeEuclideanPosition.x + " " + parentNode.nodeEuclideanPosition.y + " " + parentNode.nodeEuclideanPosition.z);

            foreach (Node child in parentNode.nodeChildren)
            {
                Vector3 childCoord = Point4d.projectAndGetVect3(child.nodeRelativeHyperbolicProjectionPosition);
                //Matrix4x4 scaler = Matrix4x4.Scale(new Vector3(globalScaler, globalScaler, globalScaler));
                //scaler.MultiplyPoint3x4(childCoord);
                //Debug.Log(rotationOfCoordinate.matrix4d);
                //rotationOfCoordinate.matrix4d.MultiplyPoint3x4(childCoord);
                //child.nodeHemsphereTheta / Math.PI * 180, 0.0, child.nodeHemspherePhi / Math.PI * 180
                //Quaternion rotation = Quaternion.Euler(parentNode.nodeHemsphereTheta / Math.PI * 180, parentNode.nodeHemspherePhi / Math.PI * 180, 0.0);
                //Quaternion rotation = Quaternion.Euler(child.nodeHemsphereTheta + parentTheta, 0.0, child.nodeHemspherePhi + parentPhi);
                //Quaternion rotation = Quaternion.Euler(parentPhi / Math.PI * 180, 0.0, parentTheta / Math.PI * 180);
                double xAngRot;
                double yAngRot;
                if (parentNode.nodeParent == null)
                {
                    xAngRot = Point4d.getRotAngleX(parentNode.nodeEuclideanPosition);
                    yAngRot = Point4d.getRotAngleY(parentNode.nodeEuclideanPosition);
                }
                else
                {
                    xAngRot = - Point4d.getRotAngleX(parentNode.nodeEuclideanPosition - parentNode.nodeParent.nodeEuclideanPosition);
                    yAngRot = Point4d.getRotAngleY(parentNode.nodeEuclideanPosition - parentNode.nodeParent.nodeEuclideanPosition);
                }

                //xAngRot += Point4d.getRotAngleX(new Point4d(childCoord.x, childCoord.y, childCoord.z));
                //yAngRot += Point4d.getRotAngleY(new Point4d(childCoord.x, childCoord.y, childCoord.z));

                //Quaternion rotation = Quaternion.Euler(xAngRot / Math.PI * 180, yAngRot / Math.PI * 180, 0.0);
                //Quaternion rotation = Quaternion.Euler(child.nodeHemspherePhi / Math.PI * 180, child.nodeHemsphereTheta / Math.PI * 180, 0.0);
                Quaternion rotation = Quaternion.Euler(0.0, 0.0, 0.0);
                //Matrix4x4 rotationMatrix = Matrix4x4.Rotate(rotation);
                //childCoord = rotationMatrix.MultiplyPoint(childCoord);

                //rotation = Quaternion.Euler(0.0, 0.0, 0.0);

                Matrix4x4 m = Matrix4x4.TRS(new Vector3(HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius) + parentNode.nodeEuclideanPosition.x, parentNode.nodeEuclideanPosition.y, parentNode.nodeEuclideanPosition.z), rotation, new Vector3(globalScaler, globalScaler, globalScaler));
                //Matrix4x4 m = Matrix4x4.TRS(new Vector3(HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius) + parentNode.nodeEuclideanPosition.x, parentNode.nodeEuclideanPosition.y, parentNode.nodeEuclideanPosition.z), rotation, new Vector3(globalScaler, globalScaler, globalScaler));

                childCoord = m.MultiplyPoint(childCoord);

                //rotation = Quaternion.Euler(parentTheta, 0.0, parentPhi);
                Debug.Log(parentNode.nodeName + " " + child.nodeName + " " + (parentTheta / Math.PI * 180) + " " + (parentPhi / Math.PI * 180));
                //m = Matrix4x4.Rotate(rotation);
                //childCoord = m.MultiplyPoint3x4(childCoord);
                


                child.nodeEuclideanPosition = new Point4d(childCoord.x, childCoord.y, childCoord.z);
                ComputeCoordinatesEuclidean(child);
            }
        }
       */
        //ComputeCoordinatesEuclideanTest(parentNode, 0.0, 0.0);
        Matrix4d I = new Matrix4d();
        I.setIdentity();
        ComputeCoordinatesEuclideanTest2(parentNode, I);
    }

    public void ComputeCoordinatesEuclideanTest(Node parentNode, double PhiCumulative, double ThetaCumulative)
    {
        /*
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            
            //Matrix4d rotationOfCoordinate = HyperbolicTransformation.buildCanonicalOrientationEuclidean(HyperbolicTransformation.ORIGIN4, parentNode.nodeEuclideanPosition);
            Point4d tempTransPoint = parentNode.nodeRelativeHyperbolicProjectionPosition;
            tempTransPoint.x += parentNode.nodeHemsphereRadius;
            Matrix4d translationMatrix = HyperbolicTransformation.buildTranslation(HyperbolicTransformation.ORIGIN4, tempTransPoint);
            Matrix4d rotationXMatrix = HyperbolicTransformation.buildXRotation(ThetaCumulative);
            Matrix4d rotationZMatrix = HyperbolicTransformation.buildZRotation(PhiCumulative);
            foreach (Node child in parentNode.nodeChildren)
            {
                Vector3 tempVect3 = translationMatrix.matrix4d.MultiplyPoint(child.nodeRelativeHyperbolicProjectionPosition.GetVector4());
                tempVect3 = rotationXMatrix.matrix4d.MultiplyPoint(tempVect3);
                tempVect3 = rotationZMatrix.matrix4d.MultiplyPoint(tempVect3);

                tempVect3 = HyperbolicTransformation.buildXRotation(child.nodeHemsphereTheta).matrix4d.MultiplyPoint(tempVect3);
                tempVect3 = HyperbolicTransformation.buildZRotation(child.nodeHemspherePhi).matrix4d.MultiplyPoint(tempVect3);
                //TODO: trans to actual euclidean position.
                Matrix4x4 scalar = Matrix4x4.Scale(new Vector3(globalScaler, globalScaler, globalScaler));
                tempVect3 = scalar.MultiplyPoint(tempVect3);
                child.nodeEuclideanPosition = new Point4d(tempVect3.x, tempVect3.y, tempVect3.z);
                
                ComputeCoordinatesEuclideanTest(child, PhiCumulative + child.nodeHemspherePhi, ThetaCumulative + child.nodeHemsphereTheta);
            }
            

        }*/
    }

    public void ComputeCoordinatesEuclideanTest2(Node parentNode, Matrix4d parentTransformation)
    {
    
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            Point4d origin = new Point4d(0, 0, 0, 1);
            foreach(Node child in parentNode.nodeChildren)
            {
                /*
                Matrix4d transformMat = HyperbolicMath.getTranslationMatrix(origin, new Point4d(parentNode.nodeHemsphereRadius, 0, 0, 1));

                transformMat = HyperbolicMath.getRotationMatrix(origin, new Point4d(0.5, 0, 0, 1), child.nodeHemsphereTheta) * transformMat;
                transformMat = HyperbolicMath.getRotationMatrix(origin, new Point4d(0, 0, 0.5, 1), child.nodeHemspherePhi) * transformMat;

                transformMat = transformMat * parentTransformation;
                transformMat.transform(child.nodeEuclideanPositionUnscaled);

                ComputeCoordinatesEuclideanTest2(child, transformMat);
                */

                Point4d childCoord = child.nodeRelativeHyperbolicProjectionPosition;
                parentTransformation.transform(childCoord);
                child.nodeEuclideanPositionUnscaled = childCoord;

                Matrix4d nextTransform = new Matrix4d();
                nextTransform = HyperbolicMath.getTranslationMatrix(origin, childCoord);

                //Vector3 childCoord = Point4d.projectAndGetVect3(child.nodeRelativeHyperbolicProjectionPosition);
                //Matrix4x4 nextRotation = new Matrix4x4();
                //Quaternion q = Quaternion.Euler(- child.nodeHemspherePhi / Math.PI * 180, child.nodeHemsphereTheta / Math.PI * 180, 0.0);
                //nextRotation = Matrix4x4.Rotate(q);
                //nextRotation *= parentTransformation;

                //Quaternion rotation = Quaternion.Euler(0.0, 0.0, 0.0);
                //Matrix4x4 m = Matrix4x4.TRS(new Vector3(HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius) + parentNode.nodeEuclideanPosition.x, parentNode.nodeEuclideanPosition.y, parentNode.nodeEuclideanPosition.z), rotation, new Vector3(globalScaler, globalScaler, globalScaler));

                //childCoord = parentTransformation.MultiplyPoint(childCoord);
                //childCoord = m.MultiplyPoint(childCoord);

                //child.nodeEuclideanPosition = new Point4d(childCoord.x, childCoord.y, childCoord.z);
                ComputeCoordinatesEuclideanTest2(child, nextTransform);
            }
        }
        
    }

    public void ComputeCoordinatesEuclideanTest2_TranslationOnly(Node parentNode, Matrix4d parentTransformation)
    {

        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            Point4d origin = new Point4d(0, 0, 0, 1);
            foreach (Node child in parentNode.nodeChildren)
            {
                Point4d childCoord = child.nodeRelativeHyperbolicProjectionPosition;
                parentTransformation.transform(childCoord);
                child.nodeEuclideanPositionUnscaled = childCoord;

                Matrix4d nextTransform = new Matrix4d();
                nextTransform = HyperbolicMath.getTranslationMatrix(origin, childCoord);

                //Vector3 childCoord = Point4d.projectAndGetVect3(child.nodeRelativeHyperbolicProjectionPosition);
                //Matrix4x4 nextRotation = new Matrix4x4();
                //Quaternion q = Quaternion.Euler(- child.nodeHemspherePhi / Math.PI * 180, child.nodeHemsphereTheta / Math.PI * 180, 0.0);
                //nextRotation = Matrix4x4.Rotate(q);
                //nextRotation *= parentTransformation;

                //Quaternion rotation = Quaternion.Euler(0.0, 0.0, 0.0);
                //Matrix4x4 m = Matrix4x4.TRS(new Vector3(HyperbolicMath.euclideanDistance(parentNode.nodeHemsphereRadius) + parentNode.nodeEuclideanPosition.x, parentNode.nodeEuclideanPosition.y, parentNode.nodeEuclideanPosition.z), rotation, new Vector3(globalScaler, globalScaler, globalScaler));

                //childCoord = parentTransformation.MultiplyPoint(childCoord);
                //childCoord = m.MultiplyPoint(childCoord);

                //child.nodeEuclideanPosition = new Point4d(childCoord.x, childCoord.y, childCoord.z);
                ComputeCoordinatesEuclideanTest2(child, nextTransform);
            }
        }

    }

    public void Scaling(Node parentNode)
    {
    /*
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }
        else
        {
            double scaler = 100f;
            foreach (Node child in parentNode.nodeChildren)
            {
                //Quaternion rotation = Quaternion.Euler(child.nodeHemsphereTheta, 0.0, child.nodeHemspherePhi);
                // Set the translation, rotation and scale parameters.
                Matrix4x4 m = Matrix4x4.Scale(new Vector3(scaler, scaler, scaler));
                Vector3 tempVect3 = m.MultiplyPoint3x4(Point4d.projectAndGetVect3(child.nodeEuclideanPosition));
                child.nodeEuclideanPosition = new Point4d(tempVect3.x, tempVect3.y, tempVect3.z);
                Scaling(child);
            }
        }
        */
    }

    public void DrawNodes()
    {
        nodesPrimitives.Add(GameObject.CreatePrimitive(PrimitiveType.Sphere));
        //nodesPrimitives[nodesPrimitives.Count - 1].transform.position = new Vector3(graphRoot.nodeEuclideanPosition.x, graphRoot.nodeEuclideanPosition.y, graphRoot.nodeEuclideanPosition.z);
        nodesPrimitives[nodesPrimitives.Count - 1].name = graphRoot.nodeName;
        DrawNodesRecursive(graphRoot);
    }

    public void DrawNodesRecursive(Node node)
    {
        foreach (Node child in node.nodeChildren)
        {
            nodesPrimitives.Add(GameObject.CreatePrimitive(PrimitiveType.Sphere));
            nodesPrimitives[nodesPrimitives.Count - 1].transform.position = new Vector3((float)child.nodeEuclideanPositionUnscaled.x * globalScaler,
                                                                                        (float)child.nodeEuclideanPositionUnscaled.y * globalScaler,
                                                                                        (float)child.nodeEuclideanPositionUnscaled.z * globalScaler);
            nodesPrimitives[nodesPrimitives.Count - 1].name = child.nodeId.ToString();
            //nodesPrimitives[nodesPrimitives.Count - 1].AddComponent<>();
            //nodesPrimitives[nodesPrimitives.Count - 1].tag = child.nodeId.ToString();
            DrawNodesRecursive(child);
        }
    }
    
    public void DrawEdges(Node parentNode)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }

        foreach(Node child in parentNode.nodeChildren)
        {
            edgeHolders.Add(new GameObject());
            edgeHolders[edgeHolders.Count - 1].AddComponent<LineRenderer>();
            LineRenderer lr = edgeHolders[edgeHolders.Count - 1].GetComponent<LineRenderer>();
            lr.material = lineMaterial; //lineMaterial
            lr.SetColors(color, color);
            lr.SetWidth(lineWidth, lineWidth);
            lr.SetPosition(0, new Vector3((float)parentNode.nodeEuclideanPositionUnscaled.x * globalScaler,
                                          (float)parentNode.nodeEuclideanPositionUnscaled.y * globalScaler,
                                          (float)parentNode.nodeEuclideanPositionUnscaled.z * globalScaler));
            lr.SetPosition(1, new Vector3((float)child.nodeEuclideanPositionUnscaled.x * globalScaler,
                                          (float)child.nodeEuclideanPositionUnscaled.y * globalScaler,
                                          (float)child.nodeEuclideanPositionUnscaled.z * globalScaler));

            DrawEdges(child);
        }
    }

    public void updateTranslation(String id)
    {
        int NodeId = Int32.Parse(id);
        Node selectedNode = graphNodesList[NodeId];
        //Node selectedNode = node;
        Matrix4d transformMat = HyperbolicMath.getTranslationMatrix(new Point4d(), selectedNode.nodeEuclideanPositionUnscaled);
        transformMat.transform(graphRoot.nodeEuclideanPositionUnscaled);
        translateAllPointByMatrix(graphRoot, transformMat);
        DrawNodes();
        DrawEdges(graphRoot);
    }

    public void translateAllPointByMatrix(Node parentNode, Matrix4d transformMat)
    {
        if (parentNode.nodeNumDecendents == 0)
        {
            return;
        }

        foreach (Node child in parentNode.nodeChildren)
        {
            transformMat.transform(child.nodeEuclideanPositionUnscaled);
            translateAllPointByMatrix(child, transformMat);
        }
    }

    public void destroyAll()
    {
        foreach(GameObject gameObject in nodesPrimitives) {
            Destroy(gameObject);
        }

        foreach(GameObject gameObject in edgeHolders)
        {
            Destroy(gameObject);
        }
    }
}


public class Node
{
    public int nodeId;
    public string nodeName;
    public int nodeLevel;
    public Point4d nodeEuclideanPosition;
    public Point4d nodeEuclideanPositionUnscaled;
    public Point4d nodeRelativeHyperbolicProjectionPosition;
    public double nodeAnglePhi;
    public double nodeAngleTheta;
    public double nodeR;
    public double nodeHemspherePhi;
    public double nodeHemsphereTheta;
    public double nodeHemsphereRadius;
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
        this.nodeEuclideanPositionUnscaled = new Point4d(0, 0, 0, 1);
        this.nodeRelativeHyperbolicProjectionPosition = new Point4d(0, 0, 0, 1);
        this.nodeAnglePhi = 0.0;
        this.nodeAngleTheta = 0.0;
        this.nodeR = 0.0;
        this.nodeHemspherePhi = 0.0;
        this.nodeHemsphereTheta = 0.0;
        this.nodeHemsphereRadius = 0.0;
        this.nodeChildren = new List<Node>();
        this.nodeParent = null;
        this.nodeNumDecendents = 0;
        this.nodeIsRoot = nodeIsRoot;
    }

    public void SetNodeEuclideanPosition(double x, double y, double z, double w)
    {
        this.nodeEuclideanPosition = new Point4d(x, y, z, w);
    }

    public void SetNodeRelativeHyperbolicProjectionPosition(double x, double y, double z, double w)
    {
        this.nodeRelativeHyperbolicProjectionPosition = new Point4d(x, y, z, w);
    }

    public void CalculateAndSetNodeEuclideanPosition()
    {
        double x = this.nodeR * Math.Sin(this.nodeAnglePhi) * Math.Cos(this.nodeAngleTheta);
        double y = this.nodeR * Math.Sin(this.nodeAnglePhi) * Math.Sin(this.nodeAngleTheta);
        double z = this.nodeR * Math.Cos(this.nodeAnglePhi);

        this.SetNodeEuclideanPosition(x, y, z, 1);
    }

    public void CalculateAndSetNodeRelativeHyperbolicProjectionPosition(double radius)
    {
        /*
        double x = radius * Math.Sin(this.nodeHemspherePhi) * Math.Cos(this.nodeHemsphereTheta);
        double y = radius * Math.Sin(this.nodeHemspherePhi) * Math.Sin(this.nodeHemsphereTheta);
        double z = radius * Math.Cos(this.nodeHemspherePhi);

        this.SetNodeRelativeHyperbolicProjectionPosition(x, y, z, 1);
        */
        
        Point4d coord = new Point4d(0, 0, 0, 1);
        Point4d origin = new Point4d(0, 0, 0, 1);
        Matrix4d transformMat = HyperbolicMath.getTranslationMatrix(origin, new Point4d(radius, 0, 0, 1));
        Matrix4d rotationMat = HyperbolicMath.getRotationMatrix(origin, new Point4d(0.5, 0, 0, 1), this.nodeHemsphereTheta);
        rotationMat *= HyperbolicMath.getRotationMatrix(coord, new Point4d(0, 0, 0.5, 1), this.nodeHemspherePhi);
        transformMat.transform(coord);
        rotationMat.transform(coord);
        coord.normalizeHomoCoord();
        this.SetNodeRelativeHyperbolicProjectionPosition(coord.x, coord.y, coord.z, coord.w);
        
    }

    public void SetHemspheres(double Phi, double Theta)
    {
        this.nodeHemspherePhi = Phi;
        this.nodeHemsphereTheta = Theta;
    }

    public void SetPolarCoordinates(double R, double Phi, double Theta)
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
    public double m00; public double m01; public double m02; public double m03;
    public double m10; public double m11; public double m12; public double m13;
    public double m20; public double m21; public double m22; public double m23;
    public double m30; public double m31; public double m32; public double m33;

    public Matrix4d(double m00 = 0.0, double m01 = 0.0, double m02 = 0.0, double m03 = 0.0,
              double m10 = 0.0, double m11 = 0.0, double m12 = 0.0, double m13 = 0.0,
              double m20 = 0.0, double m21 = 0.0, double m22 = 0.0, double m23 = 0.0,
              double m30 = 0.0, double m31 = 0.0, double m32 = 0.0, double m33 = 0.0)
    {
        this.m00 = m00; this.m01 = m01; this.m02 = m02; this.m03 = m03;
        this.m10 = m10; this.m11 = m11; this.m12 = m12; this.m13 = m13;
        this.m20 = m20; this.m21 = m21; this.m22 = m22; this.m23 = m23;
        this.m30 = m30; this.m31 = m31; this.m32 = m32; this.m33 = m33;
    }

    public void SetColumn(int index, Point4d ColumnValue)
    {
        if (index == 0)
        {
            this.m00 = ColumnValue.x;
            this.m10 = ColumnValue.y;
            this.m20 = ColumnValue.z;
            this.m30 = ColumnValue.w;
        } 
        else if (index == 1)
        {
            this.m01 = ColumnValue.x;
            this.m11 = ColumnValue.y;
            this.m21 = ColumnValue.z;
            this.m31 = ColumnValue.w;
        }
        else if (index == 2)
        {
            this.m02 = ColumnValue.x;
            this.m12 = ColumnValue.y;
            this.m22 = ColumnValue.z;
            this.m32 = ColumnValue.w;
        }
        else if (index == 3)
        {
            this.m03 = ColumnValue.x;
            this.m13 = ColumnValue.y;
            this.m23 = ColumnValue.z;
            this.m33 = ColumnValue.w;
        }
    }

    public void SetRow(int index, Point4d ColumnValue)
    {
        if (index == 0)
        {
            this.m00 = ColumnValue.x;
            this.m01 = ColumnValue.y;
            this.m02 = ColumnValue.z;
            this.m03 = ColumnValue.w;
        }
        else if (index == 1)
        {
            this.m10 = ColumnValue.x;
            this.m11 = ColumnValue.y;
            this.m12 = ColumnValue.z;
            this.m13 = ColumnValue.w;
        }
        else if (index == 2)
        {
            this.m20 = ColumnValue.x;
            this.m21 = ColumnValue.y;
            this.m22 = ColumnValue.z;
            this.m23 = ColumnValue.w;
        }
        else if (index == 3)
        {
            this.m30 = ColumnValue.x;
            this.m31 = ColumnValue.y;
            this.m32 = ColumnValue.z;
            this.m33 = ColumnValue.w;
        }
    }

    public Point4d GetRow(int index)
    {
        if (index == 0)
        {
            return new Point4d(this.m00, this.m01, this.m02, this.m03);
        }
        else if (index == 1)
        {
            return new Point4d(this.m10, this.m11, this.m12, this.m13);
        }
        else if (index == 2)
        {
            return new Point4d(this.m20, this.m21, this.m22, this.m23);
        }
        else if (index == 3)
        {
            return new Point4d(this.m30, this.m31, this.m32, this.m33);
        }

        return new Point4d();

    }

    public Point4d GetColumn(int index)
    {
        if (index == 0)
        {
            return new Point4d(this.m00, this.m10, this.m20, this.m30);
        }
        else if (index == 1)
        {
            return new Point4d(this.m01, this.m11, this.m21, this.m31);
        }
        else if (index == 2)
        {
            return new Point4d(this.m02, this.m12, this.m22, this.m32);
        }
        else if (index == 3)
        {
            return new Point4d(this.m03, this.m13, this.m23, this.m33);
        }
       
        return new Point4d();

    }

    public void transform(Point4d v)
    {
        double x = this.m00 * v.x + (this.m01 * v.y) + (this.m02 * v.z) + (this.m03 * v.w);
        double y = this.m10 * v.x + (this.m11 * v.y) + (this.m12 * v.z) + (this.m13 * v.w);
        double z = this.m20 * v.x + (this.m21 * v.y) + (this.m22 * v.z) + (this.m23 * v.w);
        double w = this.m30 * v.x + (this.m31 * v.y) + (this.m32 * v.z) + (this.m33 * v.w);

        v.x = x; v.y = y; v.z = z; v.w = w;
        v.normalizeHomoCoord();
    }



    public void rotX(double theta)
    {
        double cos_theta = Math.Cos(theta);
        double sin_theta = Math.Sin(theta);
        double neg_sin_theta = -sin_theta;

        this.m00 = 1; this.m01 = 0; this.m02 = 0; this.m03 = 0;
        this.m10 = 0; this.m11 = cos_theta; this.m12 = neg_sin_theta; this.m13 = 0;
        this.m20 = 0; this.m21 = sin_theta; this.m22 = cos_theta; this.m23 = 0;
        this.m30 = 0; this.m31 = 0; this.m32 = 0; this.m33 = 1;
    }

    public void rotY(double theta)
    {
        double cos_theta = Math.Cos(theta);
        double sin_theta = Math.Sin(theta);
        double neg_sin_theta = -sin_theta;

        this.m00 = cos_theta; this.m01 = 0; this.m02 = sin_theta; this.m03 = 0;
        this.m10 = 0; this.m11 = 1; this.m12 = 0; this.m13 = 0;
        this.m20 = neg_sin_theta; this.m21 = 0; this.m22 = cos_theta; this.m23 = 0;
        this.m30 = 0; this.m31 = 0; this.m32 = 0; this.m33 = 1;
    }

    public void rotZ(double theta)
    {
        double cos_theta = Math.Cos(theta);
        double sin_theta = Math.Sin(theta);
        double neg_sin_theta = -sin_theta;

        this.m00 = cos_theta; this.m01 = neg_sin_theta; this.m02 = 0; this.m03 = 0;
        this.m10 = sin_theta; this.m11 = cos_theta; this.m12 = 0; this.m13 = 0;
        this.m20 = 0; this.m21 = 0; this.m22 = 1; this.m23 = 0;
        this.m30 = 0; this.m31 = 0; this.m32 = 0; this.m33 = 1;
    }

    public void setIdentity()
    {
        this.m00 = 1; this.m01 = 0; this.m02 = 0; this.m03 = 0;
        this.m10 = 0; this.m11 = 1; this.m12 = 0; this.m13 = 0;
        this.m20 = 0; this.m21 = 0; this.m22 = 1; this.m23 = 0;
        this.m30 = 0; this.m31 = 0; this.m32 = 0; this.m33 = 1;
    }

    public void setHyperbolicIdentity()
    {
        this.m00 = 1; this.m01 = 0; this.m02 = 0; this.m03 = 0;
        this.m10 = 0; this.m11 = 1; this.m12 = 0; this.m13 = 0;
        this.m20 = 0; this.m21 = 0; this.m22 = 1; this.m23 = 0;
        this.m30 = 0; this.m31 = 0; this.m32 = 0; this.m33 = -1;
    }

    public static Matrix4d operator *(Matrix4d lhs, Matrix4d rhs)
    {
        Matrix4d temp = new Matrix4d();

        temp.m00 = lhs.m00 * rhs.m00 + lhs.m01 * rhs.m10 + lhs.m02 * rhs.m20 + lhs.m03 * rhs.m30;
        temp.m10 = lhs.m10 * rhs.m00 + lhs.m11 * rhs.m10 + lhs.m12 * rhs.m20 + lhs.m13 * rhs.m30;
        temp.m20 = lhs.m20 * rhs.m00 + lhs.m21 * rhs.m10 + lhs.m22 * rhs.m20 + lhs.m23 * rhs.m30;
        temp.m30 = lhs.m30 * rhs.m00 + lhs.m31 * rhs.m10 + lhs.m32 * rhs.m20 + lhs.m33 * rhs.m30;

        temp.m01 = lhs.m00 * rhs.m01 + lhs.m01 * rhs.m11 + lhs.m02 * rhs.m21 + lhs.m03 * rhs.m31;
        temp.m11 = lhs.m10 * rhs.m01 + lhs.m11 * rhs.m11 + lhs.m12 * rhs.m21 + lhs.m13 * rhs.m31;
        temp.m21 = lhs.m20 * rhs.m01 + lhs.m21 * rhs.m11 + lhs.m22 * rhs.m21 + lhs.m23 * rhs.m31;
        temp.m31 = lhs.m30 * rhs.m01 + lhs.m31 * rhs.m11 + lhs.m32 * rhs.m21 + lhs.m33 * rhs.m31;

        temp.m02 = lhs.m00 * rhs.m02 + lhs.m01 * rhs.m12 + lhs.m02 * rhs.m22 + lhs.m03 * rhs.m32;
        temp.m12 = lhs.m10 * rhs.m02 + lhs.m11 * rhs.m12 + lhs.m12 * rhs.m22 + lhs.m13 * rhs.m32;
        temp.m22 = lhs.m20 * rhs.m02 + lhs.m21 * rhs.m12 + lhs.m22 * rhs.m22 + lhs.m23 * rhs.m32;
        temp.m32 = lhs.m30 * rhs.m02 + lhs.m31 * rhs.m12 + lhs.m32 * rhs.m22 + lhs.m33 * rhs.m32;

        temp.m03 = lhs.m00 * rhs.m03 + lhs.m01 * rhs.m13 + lhs.m02 * rhs.m23 + lhs.m03 * rhs.m33;
        temp.m13 = lhs.m10 * rhs.m03 + lhs.m11 * rhs.m13 + lhs.m12 * rhs.m23 + lhs.m13 * rhs.m33;
        temp.m23 = lhs.m20 * rhs.m03 + lhs.m21 * rhs.m13 + lhs.m22 * rhs.m23 + lhs.m23 * rhs.m33;
        temp.m33 = lhs.m30 * rhs.m03 + lhs.m31 * rhs.m13 + lhs.m32 * rhs.m23 + lhs.m33 * rhs.m33;

        return temp;
    }

    public static Matrix4d operator *(Matrix4d lhs, double rhs)
    {
        Matrix4d temp = new Matrix4d();
        temp.m00 = lhs.m00 * rhs; temp.m01 = lhs.m01 * rhs; temp.m02 = lhs.m02 * rhs; temp.m03 = lhs.m03 * rhs;
        temp.m10 = lhs.m10 * rhs; temp.m11 = lhs.m11 * rhs; temp.m12 = lhs.m12 * rhs; temp.m13 = lhs.m13 * rhs;
        temp.m20 = lhs.m20 * rhs; temp.m21 = lhs.m21 * rhs; temp.m22 = lhs.m22 * rhs; temp.m23 = lhs.m23 * rhs;
        temp.m30 = lhs.m30 * rhs; temp.m31 = lhs.m31 * rhs; temp.m32 = lhs.m32 * rhs; temp.m33 = lhs.m33 * rhs;
        return temp;
    }

    public static Matrix4d operator /(Matrix4d lhs, double rhs)
    {
        Matrix4d temp = new Matrix4d();
        temp.m00 = lhs.m00 / rhs; temp.m01 = lhs.m01 / rhs; temp.m02 = lhs.m02 / rhs; temp.m03 = lhs.m03 / rhs;
        temp.m10 = lhs.m10 / rhs; temp.m11 = lhs.m11 / rhs; temp.m12 = lhs.m12 / rhs; temp.m13 = lhs.m13 / rhs;
        temp.m20 = lhs.m20 / rhs; temp.m21 = lhs.m21 / rhs; temp.m22 = lhs.m22 / rhs; temp.m23 = lhs.m23 / rhs;
        temp.m30 = lhs.m30 / rhs; temp.m31 = lhs.m31 / rhs; temp.m32 = lhs.m32 / rhs; temp.m33 = lhs.m33 / rhs;
        return temp;
    }

    public static Matrix4d operator -(Matrix4d lhs, Matrix4d rhs)
    {
        Matrix4d temp = new Matrix4d();
        temp.m00 = lhs.m00 - rhs.m00; temp.m01 = lhs.m01 - rhs.m01; temp.m02 = lhs.m02 - rhs.m02; temp.m03 = lhs.m03 - rhs.m03;
        temp.m10 = lhs.m10 - rhs.m10; temp.m11 = lhs.m11 - rhs.m11; temp.m12 = lhs.m12 - rhs.m12; temp.m13 = lhs.m13 - rhs.m13;
        temp.m20 = lhs.m20 - rhs.m20; temp.m21 = lhs.m21 - rhs.m21; temp.m22 = lhs.m22 - rhs.m22; temp.m23 = lhs.m23 - rhs.m23;
        temp.m30 = lhs.m30 - rhs.m30; temp.m31 = lhs.m31 - rhs.m31; temp.m32 = lhs.m32 - rhs.m32; temp.m33 = lhs.m33 - rhs.m33;
        return temp;
    }

    public static Matrix4d inverse(Matrix4d matrix)
    {
        double A2323 = matrix.m22 * matrix.m33 - matrix.m23 * matrix.m32;
        double A1323 = matrix.m21 * matrix.m33 - matrix.m23 * matrix.m31;
        double A1223 = matrix.m21 * matrix.m32 - matrix.m22 * matrix.m31;
        double A0323 = matrix.m20 * matrix.m33 - matrix.m23 * matrix.m30;
        double A0223 = matrix.m20 * matrix.m32 - matrix.m22 * matrix.m30;
        double A0123 = matrix.m20 * matrix.m31 - matrix.m21 * matrix.m30;
        double A2313 = matrix.m12 * matrix.m33 - matrix.m13 * matrix.m32;
        double A1313 = matrix.m11 * matrix.m33 - matrix.m13 * matrix.m31;
        double A1213 = matrix.m11 * matrix.m32 - matrix.m12 * matrix.m31;
        double A2312 = matrix.m12 * matrix.m23 - matrix.m13 * matrix.m22;
        double A1312 = matrix.m11 * matrix.m23 - matrix.m13 * matrix.m21;
        double A1212 = matrix.m11 * matrix.m22 - matrix.m12 * matrix.m21;
        double A0313 = matrix.m10 * matrix.m33 - matrix.m13 * matrix.m30;
        double A0213 = matrix.m10 * matrix.m32 - matrix.m12 * matrix.m30;
        double A0312 = matrix.m10 * matrix.m23 - matrix.m13 * matrix.m20;
        double A0212 = matrix.m10 * matrix.m22 - matrix.m12 * matrix.m20;
        double A0113 = matrix.m10 * matrix.m31 - matrix.m11 * matrix.m30;
        double A0112 = matrix.m10 * matrix.m21 - matrix.m11 * matrix.m20;

        double deter = matrix.m00 * (matrix.m11 * A2323 - matrix.m12 * A1323 + matrix.m13 * A1223)
            - matrix.m01 * (matrix.m10 * A2323 - matrix.m12 * A0323 + matrix.m13 * A0223)
            + matrix.m02 * (matrix.m10 * A1323 - matrix.m11 * A0323 + matrix.m13 * A0123)
            - matrix.m03 * (matrix.m10 * A1223 - matrix.m11 * A0223 + matrix.m12 * A0123);
        deter = 1 / deter;

        return new Matrix4d(
           deter * (matrix.m11 * A2323 - matrix.m12 * A1323 + matrix.m13 * A1223),
           deter * -(matrix.m01 * A2323 - matrix.m02 * A1323 + matrix.m03 * A1223),
           deter * (matrix.m01 * A2313 - matrix.m02 * A1313 + matrix.m03 * A1213),
           deter * -(matrix.m01 * A2312 - matrix.m02 * A1312 + matrix.m03 * A1212),
           deter * -(matrix.m10 * A2323 - matrix.m12 * A0323 + matrix.m13 * A0223),
           deter * (matrix.m00 * A2323 - matrix.m02 * A0323 + matrix.m03 * A0223),
           deter * -(matrix.m00 * A2313 - matrix.m02 * A0313 + matrix.m03 * A0213),
           deter * (matrix.m00 * A2312 - matrix.m02 * A0312 + matrix.m03 * A0212),
           deter * (matrix.m10 * A1323 - matrix.m11 * A0323 + matrix.m13 * A0123),
           deter * -(matrix.m00 * A1323 - matrix.m01 * A0323 + matrix.m03 * A0123),
           deter * (matrix.m00 * A1313 - matrix.m01 * A0313 + matrix.m03 * A0113),
           deter * -(matrix.m00 * A1312 - matrix.m01 * A0312 + matrix.m03 * A0112),
           deter * -(matrix.m10 * A1223 - matrix.m11 * A0223 + matrix.m12 * A0123),
           deter * (matrix.m00 * A1223 - matrix.m01 * A0223 + matrix.m02 * A0123),
           deter * -(matrix.m00 * A1213 - matrix.m01 * A0213 + matrix.m02 * A0113),
           deter * (matrix.m00 * A1212 - matrix.m01 * A0212 + matrix.m02 * A0112)
        );
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
    public double x;
    public double y;
    public double z;
    public double w;

    public Point4d(double x = 0, double y = 0, double z = 0, double w = 1)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    public void scale(double s)
    {
        this.x = this.x * s;
        this.y = this.y * s;
        this.z = this.z * s;
        this.w = this.w * s;
    }

    // this = s*t1 + t2
    public void scaleAdd(double s, Point4d t1, Point4d t2)
    {
        double tx = t1.x * s + (t2.x);
        double ty = t1.y * s + (t2.y);
        double tz = t1.z * s + (t2.z);
        double tw = t1.w * s + (t2.w);

        this.x = tx; this.y = ty; this.z = tz; this.w = tw;
    }

    public void normalizeHomoCoord()
    {
        double tx = this.x / (this.w);
        double ty = this.y / (this.w);
        double tz = this.z / (this.w);

        this.x = tx; this.y = ty; this.z = tz; this.w = 1.0;
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

    public void set(double x, double y, double z, double w)
    {
        this.x = x; this.y = y; this.z = z; this.w = w;
    }

    // Euclidean norm of homogeneous coordinates [and equivalent to
    // Point4d.distance(new Point4d(0, 0, 0, 0))].
    public double vectorLength(Point4d p)
    {
        double x2 = this.x * this.x;
        double y2 = this.y * this.y;
        double z2 = this.z * this.z;
        double w2 = this.w * this.w;
        return x2 + (y2) + (z2) / Math.Sqrt(w2);
    }

    // The usual vector dot product.
    public double vectorDot(Point4d v1)
    {
        double tx = this.x * (v1.x);
        double ty = this.y * (v1.y);
        double tz = this.z * (v1.z);
        double tw = this.w * (v1.w);
        return tx + (ty) + (tz) + (tw);
    }

    // The usual vector dot product computed from x, y, and z only.
    public double vectorDot3(Point4d v1)
    {
        double tx = this.x * (v1.x);
        double ty = this.y * (v1.y);
        double tz = this.z * (v1.z);
        return tx + (ty) + (tz);
    }

    // Returns the Minkowski inner product of this with y.
    public double minkowski(Point4d v)
    {
        double tx = this.x * (v.x);
        double ty = this.y * (v.y);
        double tz = this.z * (v.z);
        double tw = this.w * (v.w);
        return tx + (ty) + (tz) - (tw);
    }

    public static implicit operator Point4d(Vector4 vect)
    {
        return new Point4d(vect.x, vect.y, vect.z, vect.w);
    }
    /*
    public Vector4 GetVector4()
    {
        return new Vector4(this.x, this.y, this.z, this.w);
    }
    */
    /*
    public static Vector3 projectAndGetVect3(Point4d p)
    {
        double tx = p.x / (p.w);
        double ty = p.y / (p.w);
        double tz = p.z / (p.w);
        return new Vector3(tx, ty, tz);
    }
    */
    public static Point4d operator -(Point4d lhs, Point4d rhs)
    {
        return new Point4d(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
    }

    public static double operator *(Point4d lhs, Point4d rhs)
    {
        return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
    }

    public static Point4d operator *(Point4d lhs, double constant)
    {
        return new Point4d(lhs.x * constant, lhs.y * constant, lhs.z * constant, lhs.w * constant);
    }

    public static Point4d operator *(double constant, Point4d rhs)
    {
        return new Point4d(rhs.x * constant, rhs.y * constant, rhs.z * constant, rhs.w * constant);
    }

    public static Point4d operator +(Point4d lhs, Point4d rhs)
    {
        return new Point4d(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z, lhs.w + rhs.w);
    }

    //The coordinates should be normalized before using these affine coordinates function, because these function ignore w part of point4d
    public static Point4d subAffineCoord(Point4d point1, Point4d point2)
    {
        return new Point4d(point1.x - point2.x, point1.y - point2.y, point1.z - point2.z, 1);
    }

    public static double dotAffineCoord(Point4d point1, Point4d point2)
    {
        return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z;
    }

    public static double getAffineVectorLength(Point4d point1)
    {
        return Math.Sqrt(point1.x * point1.x + point1.y * point1.y + point1.z * point1.z);
    }

    public static Point4d makeAffineUnitVector(Point4d point1)
    {
        double length = getAffineVectorLength(point1);
        return new Point4d(point1.x / length, point1.y / length, point1.z / length, 1);
    }

    public static Matrix4d multPointWithItsTranspose(Point4d point)
    {
        Matrix4d resultMat = new Matrix4d();

        resultMat.m00 = point.x * point.x;
        resultMat.m01 = point.y * point.x;
        resultMat.m02 = point.z * point.x;
        resultMat.m03 = point.w * point.x;

        resultMat.m10 = point.x * point.y;
        resultMat.m11 = point.y * point.y;
        resultMat.m12 = point.z * point.y;
        resultMat.m13 = point.w * point.y;

        resultMat.m20 = point.x * point.z;
        resultMat.m21 = point.y * point.z;
        resultMat.m22 = point.z * point.z;
        resultMat.m23 = point.w * point.z;

        resultMat.m30 = point.x * point.w;
        resultMat.m31 = point.y * point.w;
        resultMat.m32 = point.z * point.w;
        resultMat.m33 = point.w * point.w;

        return resultMat;
    }
    /*
    public static Vector3 projectAndGetVect3FromHypToEuc(Point4d p)
    {
        Vector3 temp = projectAndGetVect3(p);
        temp.x = HyperbolicMath.euclideanDistance(temp.x);
        temp.y = HyperbolicMath.euclideanDistance(temp.y);
        temp.z = HyperbolicMath.euclideanDistance(temp.z);
        return temp;
    }
    */
    public void normalizeToAffine()
    {
        this.x /= this.w;
        this.y /= this.w;
        this.z /= this.w;
        this.w = 1;
    }

    public static double getPhiByPoint(Point4d p)
    {
        if (p.z == 0.0)
        {
            return 0.0;
        }
        return Math.Atan(Math.Sqrt(p.x * p.x + p.y * p.y) / p.z);
    }

    public static double getThetaByPoint(Point4d p)
    {
        if(p.y == 0.0)
        {
            return 0.0;
        }
        return Math.Atan(p.y / p.x);
    }

    public static double getRotAngleX(Point4d p)
    {
        double xAngRot;
        if (p.z == 0.0)
        {
            xAngRot = 0.0;
        }
        else
        {
            xAngRot = Math.Atan(p.y / p.z);
        }
        return xAngRot;
    }

    public static double getRotAngleY(Point4d p)
    {
        double yAngRot;
        if (p.z == 0.0)
        {
            yAngRot = 0.0;
        }
        else
        {
            yAngRot = Math.Atan(p.x / p.z);
        }
        return yAngRot;
    }

    public void toUnitVect()
    {
        double length;
        length = Math.Pow(this.x, 2) + Math.Pow(this.y, 2) + Math.Pow(this.z, 2);
        length = Math.Sqrt(length);

        this.x /= length;
        this.y /= length;
        this.z /= length;
    }
}

class HyperbolicMath
{

    public static double MinkowskiInnerProduct(Point4d x, Point4d y)
    {
        //double product;
        Matrix4d I31 = new Matrix4d();
        I31.setHyperbolicIdentity();
        x.normalizeHomoCoord();
        y.normalizeHomoCoord();
        Point4d tempVect = new Point4d(x.x, x.y, x.z, -x.w);
        return tempVect * y;
    }

    public static Matrix4d getReflectionMatrix(Point4d point)
    {
        Matrix4d identity = new Matrix4d();
        Matrix4d hyperIdentity = new Matrix4d();
        identity.setIdentity();
        hyperIdentity.setHyperbolicIdentity();
        point.normalizeHomoCoord();

        Matrix4d reflectionMat = new Matrix4d();
        reflectionMat = Point4d.multPointWithItsTranspose(point);
        reflectionMat *= 2;
        reflectionMat *= hyperIdentity;
        reflectionMat /= MinkowskiInnerProduct(point, point);
        reflectionMat = identity - reflectionMat;

        return reflectionMat;
    }

    public static Matrix4d getTranslationMatrix(Point4d point1, Point4d point2)
    {
        Point4d midpoint = getMidpoint(point1, point2);
        Matrix4d result = getReflectionMatrix(midpoint);
        result *= getReflectionMatrix(point1);
        return result;
    }

    public static Point4d getMidpoint(Point4d point1, Point4d point2)
    {
        point1.normalizeHomoCoord();
        point2.normalizeHomoCoord();
        double coefficientA = Math.Sqrt(MinkowskiInnerProduct(point2, point2) * MinkowskiInnerProduct(point1, point2));
        double coefficientB = Math.Sqrt(MinkowskiInnerProduct(point1, point1) * MinkowskiInnerProduct(point1, point2));

        Point4d midpoint = new Point4d();
        midpoint = point1 * coefficientA;
        midpoint += (point2 * coefficientB);
        return midpoint;
    }

    public static Matrix4d getRotationMatrix(Point4d point1, Point4d point2, double angle)
    {
        point1.normalizeHomoCoord();
        point2.normalizeHomoCoord();

        Point4d origin = new Point4d(0, 0, 0, 1);
        Point4d closestPoint = getClosestPoint(point1, point2);
        Matrix4d euclideanRotationMatrix = getEuclideanRotationMatrix(point1, point2, angle);
        Matrix4d translationMatrix = getTranslationMatrix(closestPoint, origin);
        Matrix4d inverseTranslationMatrix = Matrix4d.inverse(translationMatrix);
        Matrix4d rotationMatrix = inverseTranslationMatrix * euclideanRotationMatrix * translationMatrix;
        return rotationMatrix;
    }

    public static Point4d getClosestPoint(Point4d a, Point4d b)
    {
        Point4d aSubB = a - b;
        Point4d bSubA = b - a;
        double coefficientA = (Point4d.dotAffineCoord(a, aSubB) / Point4d.dotAffineCoord(aSubB, aSubB));
        double coefficientB = (Point4d.dotAffineCoord(b, bSubA) / Point4d.dotAffineCoord(bSubA, bSubA));
        Point4d closestPoint = coefficientA * b + coefficientB * a;
        return closestPoint;
    }

    public static Matrix4d getEuclideanRotationMatrix(Point4d point1, Point4d point2, double angle)
    {
        Point4d orientation = Point4d.subAffineCoord(point1, point2);
        orientation = Point4d.makeAffineUnitVector(orientation);
        double u1 = orientation.x;
        double u2 = orientation.y;
        double u3 = orientation.z;
        double c = Math.Cos(angle);
        double s = Math.Sin(angle);
        double c1 = 1 - c;

        return new Matrix4d(u1 * u1 + c * (1 - u1 * u1), u1 * u2 * c1 - u3 * s, u1 * u3 * c1 + u2 * s, 0,
                            u1 * u2 * c1 + u3 * s, u2 * u2 + c * (1 - u2 * u2), u2 * u3 * c1 - u1 * s, 0,
                            u1 * u3 * c1 - u2 * s, u2 * u3 * c1 + u1 * s, u3 * u3 + c * (1 - u3 * u3), 0,
                            0, 0, 0, 1);

    }

    public static double euclideanDistance(double x)
    {
        double y = (double)(Math.Cosh(x / 2.0));
        return Math.Sqrt(1.0 - 1.0 / (y * y));
    }

    public static double dotProduct(Vector3 x, Vector3 y)
    {
        return x.x * y.x + x.y * y.y + x.z * y.z;
    }

    public static double dotProduct(Point4d x, Point4d y)
    {
        return (x.x * y.x + x.y * y.y + x.z * y.z) / (x.w * y.w);
    }

    public static double vectorLength(Point4d p)
    {
        double w2 = p.w * p.w;
        return Math.Sqrt((p.x * p.x + p.y * p.y + p.z * p.z) / w2);
    }

    public static void makeUnitVector(Point4d p)
    {
        double s = p.w * vectorLength(p);
        p.x /= s;
        p.y /= s;
        p.z /= s;
        p.w = 1.0;
    }

    public static double ASinh(double value)
    {
        double temp = value + Math.Sqrt(Math.Pow(value, 2.0) + 1);
        return Math.Log(temp, Math.E);
    }
}

public class HyperbolicTransformation
{
    public static Point4d ORIGIN4 = new Point4d( 0.0, 0.0, 0.0, 1.0 );
    public static Matrix4d I4 = new Matrix4d(1.0, 0.0, 0.0, 0.0,
                                 0.0, 1.0, 0.0, 0.0,
                                 0.0, 0.0, 1.0, 0.0,
                                 0.0, 0.0, 0.0, 1.0);
    /*
    public static Matrix4d buildCanonicalOrientationEuclidean(Point4d a, Point4d b)
    {
        Point4d orientation = new Point4d(b.x - a.x, b.y - a.y, b.z - a.z);
        //double r = Math.Sqrt(orientation.x * orientation.x + orientation.y * orientation.y + orientation.z * orientation.z);
        double theta;
        double phi;
        if (orientation.x == 0.0)
        {
            theta = Math.Atan(orientation.y / orientation.x);
        }
        else
        {
            theta = 0.0;
        }

        if (orientation.z == 0.0)
        {
            phi = Math.Atan(Math.Sqrt(orientation.x * orientation.x + orientation.y * orientation.y) / orientation.z);
        }
        else
        {
            phi = 0.0;
        }
        //phi = Math.Atan(Math.Sqrt(orientation.x * orientation.x + orientation.y * orientation.y) / orientation.z);

        //Matrix4d rotationMat = new Matrix4d();
        //rotationMat.rotX(theta);
        //Matrix4d rotationMatResult = new Matrix4d();
        //rotationMatResult.rotZ(phi);
        //rotationMatResult *= rotationMat;
        //return rotationMatResult;

        Matrix4x4 rotationMat = new Matrix4x4();
        Quaternion quaternion = Quaternion.Euler(phi / Math.PI * 180, theta / Math.PI * 180, 0.0);
        rotationMat = Matrix4x4.Rotate(quaternion);
        return rotationMat;
    }*/
    /*
    public static Matrix4d buildCanonicalOrientation(Point4d a, Point4d b)
    {
        // local scratch variables; will be transformed
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

        // calculate spherical coordinates (rho, phi, theta) of pb

        // Projection to affine coordinates is necessary so that we can
        // directly reference the x, y, and z components in the following
        // calculations.
        pb.project();

        double rho = HyperbolicMath.vectorLength(pb);
        double phi = Math.Acos(pb.x / rho);
        double theta = Math.Atan2(pb.z, pb.y);

        if (phi == 0.0)
        {
            // rotate line to achieve alignment on positive x-axis 
            retval *= buildXRotation(theta);
            retval *= buildZRotation(phi);
        }

        return retval;
    }
    */
    public static Matrix4d buildXRotation(double angle)
    {
        Matrix4d m = new Matrix4d();
        m.rotX(angle);
        return m;
    }

    public static Matrix4d buildYRotation(double angle)
    {
        Matrix4d m = new Matrix4d();
        m.rotY(angle);
        return m;
    }

    public static Matrix4d buildZRotation(double angle)
    {
        Matrix4d m = new Matrix4d();
        m.rotZ(angle);
        return m;
    }

    public static Matrix4d buildTranslation(Point4d source, Point4d dest)
    {
        double aa_h = HyperbolicMath.MinkowskiInnerProduct(source, source);
        double bb_h = HyperbolicMath.MinkowskiInnerProduct(dest, dest);
        double ab_h = HyperbolicMath.MinkowskiInnerProduct(source, dest);
        double sourceScale = Math.Sqrt(bb_h * ab_h);
        double destScale = Math.Sqrt(aa_h * ab_h);
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
        double xx = point.x * point.x;
        double xy = point.x * point.y;
        double xz = point.x * point.z;
        double xw = point.x * point.w;

        double yy = point.y * point.y;
        double yz = point.y * point.z;
        double yw = point.y * point.w;

        double zz = point.z * point.z;
        double zw = point.z * point.w;

        double ww = point.w * point.w;

        double pp_h = xx + yy + zz - ww;
        double temp = -2.0f / pp_h;

        Matrix4d ppTI31 = new Matrix4d();

        ppTI31.SetRow(0, new Point4d(xx * temp + 1, xy * temp, xz * temp, -xw * temp));
        ppTI31.SetRow(1, new Point4d(xy * temp, yy * temp + 1, yz * temp, -yw * temp));
        ppTI31.SetRow(2, new Point4d(xz * temp, yz * temp, zz * temp + 1, -zw * temp));
        ppTI31.SetRow(3, new Point4d(xw * temp, yw * temp, zw * temp, -ww * temp + 1));

        return ppTI31;
    }
    /*
    private static Point4d findPivotPoint(Point4d a4, Point4d b4)
    {
        Vector3 a = new Vector3();
        Vector3 b = new Vector3();

        a = Point4d.projectAndGetVect3(a4);
        b = Point4d.projectAndGetVect3(b4);

        Vector3 a_minus_b = new Vector3(a.x - b.x, a.y - b.y, a.z - b.z);

        double p = HyperbolicMath.dotProduct(a, a_minus_b);
        double q = HyperbolicMath.dotProduct(b, a_minus_b);
        double r = HyperbolicMath.dotProduct(a_minus_b, a_minus_b);

        return new Point4d(p * b.x - q * a.x,
                   p * b.y - q * a.y,
                   p * b.z - q * a.z,
                   r);
    }
    */
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
    public double value;
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
