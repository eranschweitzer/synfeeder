
class Equipment(object):
   
    def __init__(self,path=None,fid=None):
        if path is not None:
            self.file = open(path,'w')

            s = (
            "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
            "<rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cim=\"http://iec.ch/TC57/2012/CIM-schema-cim16#\" xmlns:md=\"http://iec.ch/TC57/61970-552/ModelDescription/1#\" xmlns:entsoe=\"http://entsoe.eu/Secretariat/ProfileExtension/2#\" xmlns:neplan=\"http://www.neplan.ch#\">\n"
            )
            self.file.write(s)
        elif fid is not None:
            self.file = fid

    def close(self):
        self.file.write("</rdf:RDF>")
        self.file.close()

    def load_object(self,id=None,p=None,q=None):
        s =(
        "\t<cim:EnergyConsumer rdf:ID=\"%s\">\n" 
        "\t  <cim:EnergyConsumer.p>%0.6g</cim:EnergyConsumer.p>\n"
        "\t  <cim:EnergyConsumer.q>%0.6g</cim:EnergyConsumer.q>\n"
        "\t</cim:EnergyConsumer>\n"
        %(id,p,q))
        self.file.write(s)
    
    def conductor(self,id=None,r=0,x=0,bch=0,gch=0,r0=None,x0=None,b0ch=None,g0ch=None,length=0,basekv_id=None):
        if r0 is None:
            r0 = r
        if x0 is None:
            x0 = x
        if b0ch is None:
            b0ch = bch
        if g0ch is None:
            g0ch = gch
        s=(
        "\t<cim:ACLineSegment rdf:ID=\"%s\">\n"
        "\t  <cim:ACLineSegment.r>%0.6g</cim:ACLineSegment.r>\n"
        "\t  <cim:ACLineSegment.x>%0.6g</cim:ACLineSegment.x>\n"
        "\t  <cim:ACLineSegment.bch>%0.6g</cim:ACLineSegment.bch>\n"
        "\t  <cim:ACLineSegment.gch>%0.6g</cim:ACLineSegment.gch>\n"
        "\t  <cim:Conductor.length>%0.6g</cim:Conductor.length>\n"
        "\t  <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#%s\"/>\n"
        "\t  <cim:ACLineSegment.r0>%0.6g</cim:ACLineSegment.r0>\n"
        "\t  <cim:ACLineSegment.x0>%0.6g</cim:ACLineSegment.x0>\n"
        "\t  <cim:ACLineSegment.b0ch>%0.6g</cim:ACLineSegment.b0ch>\n"
        "\t  <cim:ACLineSegment.g0ch>%0.6g</cim:ACLineSegment.g0ch>\n"
        "\t</cim:ACLineSegment>\n" 
        %(id,r,x,bch,gch,length,basekv_id,r0,x0,b0ch,g0ch)
        )
        self.file.write(s)
    
    def terminal(self,id=None,sequence_num=None,equipment_id=None):
        s = (
        "\t<cim:Terminal rdf:ID=\"%s\">\n"
        "\t  <cim:Terminal.sequenceNumber>%d</cim:Terminal.sequenceNumber>\n"
        "\t  <cim:Terminal.ConductingEquipment rdf:resource=\"#%s\"/>\n"
        "\t</cim:Terminal>\n"
        %(id,sequence_num,equipment_id)
        )
        self.file.write(s)
    
    def external_network(self,id=None,sub_id=None):
        s = (
        "\t<cim:ExternalNetworkInjection rdf:ID=\"%s\">\n"
        "\t  <cim:Equipment.EquipmentContainer rdf:resource=\"#%s\"/>\n"
        "\t  <cim:ExternalNetworkInjection.voltageFactor>1.000000</cim:ExternalNetworkInjection.voltageFactor>\n"
        "\t  <cim:ExternalNetworkInjection.referencePriority>1</cim:ExternalNetworkInjection.referencePriority>\n"
        "\t</cim:ExternalNetworkInjection>\n"
        %(id,sub_id)
        )
        self.file.write(s)

    def base_voltage(self,unoms=None):
        """ create a map between voltage levels and integer ids 
        Then create code for the base voltages"""
        self.kvid = dict( zip(unoms, range(len(unoms)) ) )
        for kv, i in self.kvid.items():
            s = (
            "\t<cim:BaseVoltage rdf:ID=\"%s\">\n"
            "\t  <cim:BaseVoltage.nominalVoltage>%0.2f</cim:BaseVoltage.nominalVoltage>\n"
            "\t</cim:BaseVoltage>\n"
            %('_basekv_class_%d' %(i),kv)
            )
            self.file.write(s)

    def power_transformer_end(self,id=None,snom=None,unom=None,r=1e-6,x=1e-6,b=0,g=0,r0=1e-6,x0=1e-6,b0=0,g0=0,terminal=None,clocknum=0,winding=None):
        s = (
        "\t<cim:PowerTransformerEnd rdf:ID=\"%s\">\n"
        "\t  <cim:PowerTransformerEnd.ratedS>%0.6g</cim:PowerTransformerEnd.ratedS>\n"
        "\t  <cim:PowerTransformerEnd.ratedU>%0.6g</cim:PowerTransformerEnd.ratedU>\n"
        "\t  <cim:PowerTransformerEnd.r>%0.6g</cim:PowerTransformerEnd.r>\n"
        "\t  <cim:PowerTransformerEnd.x>%0.6g</cim:PowerTransformerEnd.x>\n"
        "\t  <cim:PowerTransformerEnd.b>%0.6g</cim:PowerTransformerEnd.b>\n"
        "\t  <cim:PowerTransformerEnd.g>%0.6g</cim:PowerTransformerEnd.g>\n"
        "\t  <cim:PowerTransformerEnd.r0>%0.6g</cim:PowerTransformerEnd.r0>\n"
        "\t  <cim:PowerTransformerEnd.x0>%0.6g</cim:PowerTransformerEnd.x0>\n"
        "\t  <cim:PowerTransformerEnd.b0>%0.6g</cim:PowerTransformerEnd.b0>\n"
        "\t  <cim:PowerTransformerEnd.g0>%0.6g</cim:PowerTransformerEnd.g0>\n"
		"\t  <cim:PowerTransformerEnd.phaseAngleClock>%d</cim:PowerTransformerEnd.phaseAngleClock>\n"
		"\t  <cim:PowerTransformerEnd.connectionKind rdf:resource=\"http://iec.ch/TC57/2013/CIM-schema-cim16#WindingConnection.%s\"/>\n"
        "\t  <cim:TransformerEnd.Terminal rdf:resource=\"#%s\"/>\n"
        "\t</cim:PowerTransformerEnd>\n"
        %(id,snom,unom,r,x,b,g,r0,x0,b0,g0,clocknum,winding,terminal)
        )
        self.file.write(s)

    def power_transformer(self,id=None,wf_id=None,wt_id=None):
        s = (
        "\t<cim:PowerTransformer rdf:ID=\"%s\">\n"
        "\t  <cim:PowerTransformer.PowerTransformerEnd rdf:resource=\"#%s\"/>\n"
        "\t  <cim:PowerTransformer.PowerTransformerEnd rdf:resource=\"#%s\"/>\n"
        "\t</cim:PowerTransformer>\n"
        %(id,wf_id,wt_id)
        )
        self.file.write(s)

    def operational_limit_set(self,id=None, equipment_id="", terminal_id=""):
        s = (
        "\t<cim:OperationalLimitSet rdf:ID=\"%s\">\n"
        "\t  <cim:OperationalLimitSet.Equipment rdf:resource=\"#%s\"/>\n"
        "\t  <cim:OperationalLimitSet.Terminal rdf:resource=\"#%s\"/>\n"
        "\t</cim:OperationalLimitSet>\n"
        %(id,equipment_id,terminal_id)
        )
        self.file.write(s)
    
    def current_limit(self, id=None, val=0, set_id=None, type_id=None):
        s = (
        "\t<cim:CurrentLimit rdf:ID=\"%s\">\n"
        "\t  <cim:CurrentLimit.value>%0.3g</cim:CurrentLimit.value>\n"
        "\t  <cim:OperationalLimit.OperationalLimitSet rdf:resource=\"#%s\"/>\n"
        "\t  <cim:OperationalLimit.OperationalLimitType rdf:resource=\"#%s\"/>\n"
        "\t</cim:CurrentLimit>\n"
        %(id,val,set_id,type_id)
        )
        self.file.write(s)

    def operational_limit_type(self, id=None):
        """ for now just generates an absolute value limit type """

        s = (
        "\t<cim:OperationalLimitType rdf:ID=\"%s\">\n"
        "\t  <cim:OperationalLimitType.direction rdf:resource=\"http://iec.ch/TC57/2012/CIM-schema-cim16#OperationalLimitDirectionKind.absoluteValue\"/>\n"
        "\t</cim:OperationalLimitType>\n"
        %(id)
        )
        self.file.write(s)

class Topology(object):
    def __init__(self,path=None,fid=None):
        if path is not None:
            self.file = open(path,'w')

            s = (
            "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
            "<rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cim=\"http://iec.ch/TC57/2012/CIM-schema-cim16#\" xmlns:md=\"http://iec.ch/TC57/61970-552/ModelDescription/1#\" xmlns:entsoe=\"http://entsoe.eu/Secretariat/ProfileExtension/2#\" xmlns:neplan=\"http://www.neplan.ch#\">\n"
            )
            self.file.write(s)
        elif fid is not None:
            self.file = fid

    def close(self):
        self.file.write("</rdf:RDF>")
        self.file.close()

    def top_node(self,id=None, basekv_id=None):
        s = (
        "\t<cim:TopologicalNode rdf:ID=\"%s\">\n"
        "\t  <cim:TopologicalNode.BaseVoltage rdf:resource=\"#%s\"/>\n"
        "\t</cim:TopologicalNode>\n"
        %(id,'_basekv_class_%d' %(basekv_id))
        )
        self.file.write(s)

    def terminal_connect(self,terminal_id=None,top_node_id=None,status=None):
        s = (
        "\t<cim:Terminal rdf:about=\"#%s\">\n"
        "\t  <cim:Terminal.connected>%s</cim:Terminal.connected>\n"
        "\t  <cim:Terminal.TopologicalNode rdf:resource=\"#%s\"/>\n"
        "\t</cim:Terminal>\n"
        %(terminal_id,str(status).casefold(),top_node_id)
        )
        self.file.write(s)


